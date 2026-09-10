#!/usr/bin/env bash
#
# Stage the macrophage bulk RNA-seq FASTQs for GEO upload, on the HPC.
#
#   ./stage_geo_upload.sh
#   SRC=/other/fastqs STAGE=~/scratch/other ./stage_geo_upload.sh
#
# Builds a flat folder of symlinks, verifies it, and writes a checksum file that
# `build_geo_metadata.py --md5` can consume.  Idempotent -- re-running relinks in place.
# Runbook: rna/GEO_SUBMISSION.md, Phases 1-2.

set -euo pipefail

SRC=${SRC:-/lustre/fs4/vino_lab/store/globus/msk_jahan/Macrophage_RNAseq_globus/231027_A00815_0704_BHKCW3DMXY/fastqs}
STAGE=${STAGE:-$HOME/scratch/macrophage_geo_sub}
MD5OUT=${MD5OUT:-$HOME/scratch/geo_md5sums.txt}

EXPECT_FILES=54    # 27 samples x R1/R2
EXPECT_SAMPLES=27  # TZ10/TZ11/TZ12 x 9 conditions; TZ3 excluded

fails=0
pass() { printf '  PASS  %s\n' "$1"; }
fail() { printf '  FAIL  %s\n' "$1"; fails=$((fails + 1)); }

check_empty() {  # check_empty <description> <output>
	if [ -z "$2" ]; then pass "$1"; else
		fail "$1"
		printf '%s\n' "$2" | sed 's/^/          /'
	fi
}

# Step 0 -- refuse the wrong source directory

echo "source: $SRC"
[ -d "$SRC" ] || { echo "ERROR: not a directory: $SRC" >&2; exit 1; }

n_src=$(find "$SRC" -maxdepth 1 -name '*.fastq.gz' | wc -l | tr -d ' ')
if [ "$n_src" -ne "$EXPECT_FILES" ]; then
	cat >&2 <<-EOF
		ERROR: $SRC holds $n_src .fastq.gz files, expected $EXPECT_FILES.

		The complete delivery from the sequencing core is
		  /lustre/fs4/vino_lab/store/globus/msk_jahan/Macrophage_RNAseq_globus/231027_A00815_0704_BHKCW3DMXY/fastqs
		Do NOT stage from
		  /lustre/fs4/vino_lab/store/tzhang01/20231103_Macrophage_RNA/fastqs
		-- that is a partial copy missing TZ11_5-TZ11_9, TZ12_1 and TZ12_9.

		Nothing was staged.
	EOF
	exit 1
fi
echo "        $n_src FASTQs found -- complete delivery"

# Step 1 -- build the symlink farm
echo
echo "staging: $STAGE"
mkdir -p "$STAGE"

# Only symlinks are cleared: the processed matrices, scp'd in later, are left alone.
find "$STAGE" -maxdepth 1 -name '*.fastq.gz' -type l -exec rm -f {} +

linked=0
for f in "$SRC"/*.fastq.gz; do
	b=$(basename "$f")
	# TZ3 is excluded from the submission (only 8 of 9 conditions).
	case "$b" in TZ3_*) continue ;; esac
	ln -sfn "$f" "$STAGE/$b"
	linked=$((linked + 1))
done
echo "        $linked symlinks"

# Step 2 -- verify before anything leaves the cluster
echo
echo "checks:"

names=$(cd "$STAGE" && ls *.fastq.gz)

n_staged=$(printf '%s\n' "$names" | wc -l | tr -d ' ')
[ "$n_staged" -eq "$EXPECT_FILES" ] &&
	pass "$EXPECT_FILES files staged" ||
	fail "$n_staged files staged, expected $EXPECT_FILES"

# With -L, find resolves working links (type f) and leaves dangling ones as type l.
check_empty "no dangling symlinks" "$(find -L "$STAGE" -maxdepth 1 -type l)"

n_samples=$(printf '%s\n' "$names" | sed 's/_S[0-9]*_R[12]_001\.fastq\.gz$//' | sort -u | wc -l | tr -d ' ')
[ "$n_samples" -eq "$EXPECT_SAMPLES" ] &&
	pass "$EXPECT_SAMPLES samples" ||
	fail "$n_samples samples, expected $EXPECT_SAMPLES"

check_empty "every R1 has a matching R2" \
	"$(printf '%s\n' "$names" | sed 's/_R[12]_/_R#_/' | sort | uniq -c | awk '$1!=2')"

# S-numbers must run 1..27, twice each, matching sample_manifest.tsv.
check_empty "S1-S$EXPECT_SAMPLES each appear twice" \
	"$(printf '%s\n' "$names" |
		sed 's/.*_S\([0-9]*\)_R[12]_001\.fastq\.gz$/\1/' |
		sort -n | uniq -c |
		awk -v n="$EXPECT_SAMPLES" '$1!=2 || $2<1 || $2>n')"

check_empty "no TZ3 samples" "$(printf '%s\n' "$names" | grep '^TZ3_' || true)"

# GEO rejects anything outside [A-Za-z0-9._-] in a filename.
check_empty "filenames are GEO-safe" "$(printf '%s\n' "$names" | grep '[^A-Za-z0-9._-]' || true)"

# Step 3 -- checksums, without re-hashing a quarter-terabyte
echo
echo "checksums:"

core=""
for m in "$SRC"/*.md5 "$SRC"/../*.md5; do
	if [ -r "$m" ]; then core=$m; break; fi
done

if [ -z "$core" ]; then
	cat <<-EOF
		  SKIP  no readable .md5 alongside the FASTQs

		        The core ships 231027_A00815_0704_BHKCW3DMXY.md5; if it is there but
		        unreadable, ask the owner for read access -- that is far cheaper than
		        re-hashing.  To compute them yourself, submit this rather than running
		        it on a login node:

		          sbatch --job-name=geo_md5 --time=8:00:00 --mem=4G \\
		                 --wrap "cd $STAGE && md5sum *.fastq.gz > $MD5OUT"
	EOF
else
	echo "        $core"
	# Reduce the core's file to just the staged names, as `<md5>  <basename>`.
	printf '%s\n' "$names" >"$STAGE/.staged_names"
	awk 'NR==FNR { want[$0]; next }
	     NF==2   { d=$1; n=$2; if (length($1)!=32) { d=$2; n=$1 }
	               sub(/.*\//,"",n); if (n in want) print d"  "n }' \
		"$STAGE/.staged_names" "$core" >"$MD5OUT"
	rm -f "$STAGE/.staged_names"

	n_md5=$(wc -l <"$MD5OUT" | tr -d ' ')
	[ "$n_md5" -eq "$EXPECT_FILES" ] &&
		pass "$EXPECT_FILES checksums -> $MD5OUT" ||
		fail "$n_md5 checksums in $MD5OUT, expected $EXPECT_FILES (stale or partial .md5?)"

	# Spot-check the two smallest files rather than re-hashing all 54.
	if [ "$n_md5" -eq "$EXPECT_FILES" ]; then
		smallest=$(cd "$STAGE" && ls -LS *.fastq.gz | tail -2)
		for b in $smallest; do
			want=$(awk -v n="$b" '$2==n {print $1}' "$MD5OUT")
			got=$(md5sum "$STAGE/$b" | awk '{print $1}')
			[ "$want" = "$got" ] &&
				pass "spot-check $b" ||
				fail "spot-check $b: .md5 says $want, file hashes to $got"
		done
	fi
fi

# Step 4 -- what is still missing
cat <<-EOF

	processed data files (not on the cluster -- do these on the Mac, from the repo root):

	  python src/data_deposition/build_processed_files.py            # both matrices + all checks
	  gzip -k rna/counts_27samples.tsv rna/normalized_counts.tsv
	  scp rna/counts_27samples.tsv.gz rna/normalized_counts.tsv.gz <login-host>:$STAGE/

	then bring the raw checksums down, add the two processed ones, and fill the sheet:

	  scp <login-host>:$MD5OUT src/data_deposition/
	  (cd rna && md5 -r counts_27samples.tsv.gz normalized_counts.tsv.gz) | sed 's/ /  /' \\
	      >>src/data_deposition/$(basename "$MD5OUT")
	  cd src/data_deposition && python build_geo_metadata.py --md5 $(basename "$MD5OUT")
EOF

echo
if [ "$fails" -eq 0 ]; then
	echo "all checks passed -- $STAGE is ready to upload"
else
	echo "$fails check(s) FAILED -- do not upload until these are resolved" >&2
	exit 1
fi
