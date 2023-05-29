zgrep -P -o -m 1 '[^,]+\.raw' *.gz | perl -pe 's/^.*?([A-Za-z0-9]+)_(R\d).*:(.*)/$1\t$2\t$3/;' | sort | uniq > sample_description.csv

zcat *prio1*.gz | head -1 > prio1.csv
zcat *prio1*.gz | grep -P -v '^Fraction' >> prio1.csv

zcat *prio2*.gz | head -1 > prio2.csv
zcat *prio2*.gz | grep -P -v '^Fraction' >> prio2.csv

ls | grep -P '_R\d\.csv\..*.gz' | xargs zcat | head -1 > FrameShift.csv
ls | grep -P '_R\d\.csv\..*.gz' | xargs zcat | grep -P -v '^Fraction' >> FrameShift.csv

