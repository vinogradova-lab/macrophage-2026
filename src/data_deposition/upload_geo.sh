#!/usr/bin/env bash
#
# Upload the staged macrophage bulk RNA-seq folder to GEO's private FTP server, from the HPC.
#
#   export GEO_FOLDER='uploads/<your_email>_<random>'   # from the logged-in GEO FTP page
#   ./upload_geo.sh --smoke     # preflight + one small file: proves the path works, ~10 s
#   ./upload_geo.sh --sbatch    # submit the full transfer to Slurm -- the recommended route
#   ./upload_geo.sh --mirror    # run the transfer here and now (inside tmux or a Slurm job)
#   ./upload_geo.sh --verify    # remote listing vs local: count, sizes, no empty files
#
# TWO values are yours to supply, and neither is the server address:
#
#   1. the geoftp passcode -> ~/.netrc, and nowhere else:
#
#          cat >> ~/.netrc <<'EOF'
#          machine ftp-private.ncbi.nlm.nih.gov
#          login geoftp
#          password <passcode>
#          EOF
#          chmod 600 ~/.netrc
#
#   2. your personalized upload folder -> $GEO_FOLDER, exported before running this.
#
# Runbook: rna/GEO_SUBMISSION.md, Phases 4-5.

set -euo pipefail

HOST=ftp-private.ncbi.nlm.nih.gov
STAGE=${STAGE:-$HOME/scratch/macrophage_geo_sub}
SUB=${SUB:-macrophage_TLR_RNAseq}   # one subfolder per submission -- GEO requires a fresh one
LOG=${LOG:-$HOME/scratch/geo_upload.log}   # outside $STAGE, so it can never be swept into the upload
NETRC=$HOME/.netrc

# Slurm, for --sbatch.  node281 (partition vino_a) reaches ftp-private:21 and has lftp.
SLURM_ACCOUNT=${SLURM_ACCOUNT:-vino_condo_bank}
SLURM_PARTITION=${SLURM_PARTITION:-vino_a}
SLURM_TIME=${SLURM_TIME:-1-00:00:00}   # 24 h; vino_a allows 7 days, and --continue resumes anyway
SLURM_CPUS=${SLURM_CPUS:-2}
SLURM_MEM=${SLURM_MEM:-4G}

EXPECT_FILES=56    # 54 FASTQ symlinks + 2 processed matrices
EXPECT_FASTQ=54

fails=0
pass() { printf '  PASS  %s\n' "$1"; }
fail() { printf '  FAIL  %s\n' "$1"; fails=$((fails + 1)); }

usage() {
	# The header comment block is the help text.
	awk 'NR>2 { if (!/^#/) exit; sub(/^# ?/, ""); print }' "$0"
	exit "${1:-1}"
}

mode=${1:-}
case "$mode" in
	--smoke | --sbatch | --mirror | --verify) ;;
	-h | --help) usage 0 ;;
	*) echo "ERROR: need exactly one of --smoke, --sbatch, --mirror, --verify" >&2; echo >&2; usage ;;
esac

# lftp settings shared by every connection below.  This string MUST stay on one line:
# given an -e argument containing newlines, lftp skips the implicit `open <site>` and every
# command fails with "Not connected".  Semicolons, not newlines.
SETTINGS="set ftp:passive-mode true; set net:timeout 30; set net:max-retries 10; set net:reconnect-interval-base 15; set mirror:parallel-transfer-count 4; set cmd:fail-exit true;"

# run_lftp <command-string> -- one line in, settings and `bye` added.  $HOST is passed
# without a protocol prefix; writing ftp://$HOST would skip ~/.netrc and log in anonymously.
run_lftp() {
	lftp -e "$SETTINGS $1; bye" "$HOST"
}

# Preflight -- every mode, before a single byte moves
echo "preflight:"

# (1) The upload folder.  GEO prints it both as "/uploads/<id>_<random>" and as
#     "cd uploads/<id>_<random>", so strip a leading or trailing slash rather than reject.
if [ -n "${GEO_FOLDER:-}" ]; then
	GEO_FOLDER=${GEO_FOLDER#/}
	GEO_FOLDER=${GEO_FOLDER%/}
fi

if [ -z "${GEO_FOLDER:-}" ]; then
	fail "GEO_FOLDER is not set -- export GEO_FOLDER='uploads/<your_id>_<random>' (from the logged-in GEO FTP page)"
elif [ "${GEO_FOLDER#uploads/}" = "$GEO_FOLDER" ]; then
	fail "GEO_FOLDER='$GEO_FOLDER' does not start with 'uploads/' -- copy the path from the GEO page verbatim"
elif [ "${GEO_FOLDER%_*}" = "$GEO_FOLDER" ]; then
	fail "GEO_FOLDER='$GEO_FOLDER' has no '_<random>' suffix -- looks truncated"
else
	pass "GEO_FOLDER=$GEO_FOLDER"
fi

# (2) Credentials.  A missing .netrc stanza logs in anonymously and fails on the first write.
if [ ! -f "$NETRC" ]; then
	fail "$NETRC does not exist -- create it (see the header of this script); it is yours to make, not the cluster's"
else
	perm=$(stat -c '%a' "$NETRC")
	[ "$perm" = "600" ] &&
		pass "$NETRC is mode 600" ||
		fail "$NETRC is mode $perm -- run: chmod 600 $NETRC  (lftp ignores a group/world-readable .netrc)"

	if grep -qE "^[[:space:]]*machine[[:space:]]+$HOST([[:space:]]|$)" "$NETRC"; then
		pass "$NETRC has a '$HOST' stanza"
		if grep -qE '^[[:space:]]*password[[:space:]]+(YOUR_PASSCODE_HERE|PASTE_GEO_PASSCODE_HERE|<passcode>)' "$NETRC"; then
			fail "$NETRC still holds the placeholder passcode -- paste the real one from the GEO page"
		else
			pass "passcode line is filled in"
		fi
	else
		fail "$NETRC has no 'machine $HOST' stanza -- lftp would log in anonymously"
	fi
fi

# (3) The staged folder.  A dangling symlink uploads as nothing; a stray file gets submitted.
if [ ! -d "$STAGE" ]; then
	fail "$STAGE is not a directory -- run stage_geo_upload.sh first"
else
	n_all=$(find "$STAGE" -maxdepth 1 -mindepth 1 | wc -l | tr -d ' ')
	n_fq=$(find "$STAGE" -maxdepth 1 -name '*.fastq.gz' | wc -l | tr -d ' ')
	[ "$n_all" -eq "$EXPECT_FILES" ] &&
		pass "$EXPECT_FILES entries staged" ||
		fail "$n_all entries in $STAGE, expected $EXPECT_FILES (54 FASTQs + 2 matrices)"
	[ "$n_fq" -eq "$EXPECT_FASTQ" ] &&
		pass "$EXPECT_FASTQ FASTQs" ||
		fail "$n_fq FASTQs in $STAGE, expected $EXPECT_FASTQ"

	dangling=$(find -L "$STAGE" -maxdepth 1 -type l)
	[ -z "$dangling" ] &&
		pass "no dangling symlinks" ||
		{ fail "dangling symlinks"; printf '%s\n' "$dangling" | sed 's/^/          /'; }
fi

# (4) Egress.
if command -v nc >/dev/null && nc -z -w 10 "$HOST" 21 2>/dev/null; then
	pass "port 21 to $HOST reachable"
else
	fail "cannot reach $HOST:21 -- outbound FTP looks blocked from $(hostname -s)"
fi

if [ "$fails" -ne 0 ]; then
	echo
	echo "$fails preflight check(s) FAILED -- nothing was uploaded" >&2
	exit 1
fi

total=$(find -L "$STAGE" -maxdepth 1 -type f -printf '%s\n' | awk '{s+=$1} END {print s}')
echo
echo "target:  $HOST  ->  $GEO_FOLDER/$SUB"
echo "source:  $STAGE  ($EXPECT_FILES files, $(numfmt --to=iec "$total" 2>/dev/null || echo "$total bytes"))"
echo

# --sbatch -- hand the transfer to Slurm and stop caring about the SSH session
if [ "$mode" = "--sbatch" ]; then
	command -v sbatch >/dev/null || { echo "ERROR: no sbatch on $(hostname -s)" >&2; exit 1; }
	self=$(readlink -f "$0")

	jid=$(sbatch --parsable \
		--account="$SLURM_ACCOUNT" \
		--partition="$SLURM_PARTITION" \
		--time="$SLURM_TIME" \
		--cpus-per-task="$SLURM_CPUS" \
		--mem="$SLURM_MEM" \
		--job-name=geo_upload \
		--output="$HOME/scratch/geo_upload_%j.out" \
		--export="ALL,GEO_FOLDER=$GEO_FOLDER,SUB=$SUB,STAGE=$STAGE,LOG=$LOG" \
		--wrap "$self --mirror")

	cat <<-EOF
		submitted job $jid to $SLURM_PARTITION (account $SLURM_ACCOUNT, walltime $SLURM_TIME)

		  squeue -j $jid                              # queued / running
		  tail -f $HOME/scratch/geo_upload_$jid.out   # live progress
		  scancel $jid                                # stop it; --continue resumes on the next run

		You can log out.  When it finishes:

		  GEO_FOLDER='$GEO_FOLDER' $self --verify
	EOF
	exit 0
fi

# --smoke -- prove the path end to end before committing hours to it
if [ "$mode" = "--smoke" ]; then
	# Passive data ports are a separate firewall rule from port 21: a hang after `cd`
	# succeeds is that failure.
	echo "smoke test: creating $GEO_FOLDER/$SUB and uploading one 1.9 MB file"
	echo
	# Separate connection, without fail-exit: lftp errors on an existing dir even under -p.
	lftp -e "set ftp:passive-mode true; mkdir -p $GEO_FOLDER/$SUB; bye" "$HOST" || true

	run_lftp "cd $GEO_FOLDER/$SUB; put $STAGE/counts_27samples.tsv.gz; cls -1 -s --filesize --block-size=1 -B"
	cat <<-EOF

		smoke test passed.  The file above is real -- --mirror will find it already present
		and skip it, so nothing here is wasted work.

		Now start the full transfer on a compute node, so it survives your SSH session and
		keeps hours of network I/O off the login node:

		  GEO_FOLDER='$GEO_FOLDER' $0 --sbatch
	EOF
	exit 0
fi

# --mirror -- the full transfer
if [ "$mode" = "--mirror" ]; then
	# Refuse to run somewhere a dropped SSH session would kill it.
	if [ -z "${SLURM_JOB_ID:-}" ] && [ -z "${TMUX:-}" ] && [ "${FORCE:-0}" != "1" ]; then
		cat >&2 <<-EOF
			ERROR: this is a login-node shell outside tmux.

			The transfer moves ~389 GB over several hours; a dropped SSH session would kill it,
			and a login node is the wrong place for that much sustained network I/O.  Either:

			  GEO_FOLDER='$GEO_FOLDER' $0 --sbatch      # recommended -- runs on a compute node

			or keep it here, in a session that survives disconnection:

			  tmux new -s geo
			  GEO_FOLDER='$GEO_FOLDER' $0 --mirror      # ctrl-b d detaches, tmux attach -t geo returns

			FORCE=1 overrides this if you have arranged the session yourself (nohup, srun --pty
			inside tmux).
		EOF
		exit 1
	fi

	if [ -n "${SLURM_JOB_ID:-}" ]; then
		echo "running as Slurm job $SLURM_JOB_ID on $(hostname -s)"
	fi

	echo "mirroring -- log: $LOG"
	echo "  (ctrl-b d detaches tmux; the transfer keeps running)"
	echo
	# -R uploads, -L dereferences the symlink farm, --continue resumes partial files.
	# --no-perms is required: GEO refuses mirror's post-upload CHMOD ("550 Permission
	# denied"), which under cmd:fail-exit would fail a job whose bytes all arrived.
	run_lftp "cd $GEO_FOLDER/$SUB; mirror -R -L --continue --no-perms --verbose $STAGE ." 2>&1 | tee -a "$LOG"

	cat <<-EOF

		transfer finished.  Confirm it before submitting metadata:

		  GEO_FOLDER='$GEO_FOLDER' $0 --verify
	EOF
	exit 0
fi

# --verify -- the remote copy against the local one, file by file

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

run_lftp "cd $GEO_FOLDER/$SUB; cls -1 -s --filesize --block-size=1 -B" |
	awk 'NF==2 {print $2" "$1}' | sort >"$tmp/remote"

find -L "$STAGE" -maxdepth 1 -type f -printf '%f %s\n' | sort >"$tmp/local"

echo "checks:"

n_remote=$(wc -l <"$tmp/remote" | tr -d ' ')
[ "$n_remote" -eq "$EXPECT_FILES" ] &&
	pass "$EXPECT_FILES files present remotely" ||
	fail "$n_remote files in $GEO_FOLDER/$SUB, expected $EXPECT_FILES"

empty=$(awk '$2==0 {print $1}' "$tmp/remote")
[ -z "$empty" ] &&
	pass "no zero-byte files" ||
	{ fail "zero-byte files remotely"; printf '%s\n' "$empty" | sed 's/^/          /'; }

diffs=$(diff "$tmp/local" "$tmp/remote" || true)
if [ -z "$diffs" ]; then
	pass "every filename and size matches local, byte for byte"
else
	fail "local and remote disagree (< local, > remote)"
	printf '%s\n' "$diffs" | sed 's/^/          /'
fi

sum_l=$(awk '{s+=$2} END {print s+0}' "$tmp/local")
sum_r=$(awk '{s+=$2} END {print s+0}' "$tmp/remote")
[ "$sum_l" = "$sum_r" ] &&
	pass "total bytes match ($sum_l)" ||
	fail "total bytes differ: local $sum_l, remote $sum_r"

echo
if [ "$fails" -eq 0 ]; then
	cat <<-EOF
		all checks passed -- the upload is complete and matches $STAGE.

		The upload is NOT a submission.  Files sitting on the FTP server with no metadata are
		invisible to curators and are eventually purged.  Finish at

		  https://submit.ncbi.nlm.nih.gov/geo/submission/   -> Submit Metadata

		upload src/data_deposition/seq_template_filled.xlsx and give the subfolder name: $SUB

		Do not re-upload or change files after that -- corrections go to geo@ncbi.nlm.nih.gov.
	EOF
else
	echo "$fails check(s) FAILED -- re-run --mirror (--continue resumes) before submitting metadata" >&2
	exit 1
fi
