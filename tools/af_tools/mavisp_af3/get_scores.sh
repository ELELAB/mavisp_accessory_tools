. /usr/local/envs/py310/bin/activate

set -euo pipefail

ZIPFILE=$1

unzip $ZIPFILE

./extract_pae . -v -a AF3

for f in $(ls fold_*.cif); do
	BeEM -p=${f%%.*} $f
done

rm -f scores.txt
for f in *summary_confidences*; do
	modeln=$(echo $f | sed -E 's/.*summary_confidences_([0-9]+).json/\1/')
	echo "model $modeln" >> scores.txt
	grep \"ptm\" $f >> scores.txt
	grep \"iptm\" $f >> scores.txt
done

