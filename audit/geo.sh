for g in GSE113230 GSE58135 GSE133998 GSE77509 GSE130078 GSE248612 GSE221107 GSE229705 GSE148036 GSE121842 GSE136569 GSE143797 GSE108735 GSE171011 GSE97239 GSE192410 GSE192849; do
  s=$(curl -s -m 40 "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$g&targ=self&form=text&view=brief" | grep -m1 "Series_sample_id" | sed 's/.*= //' | tr -d '\r ')
  sleep 2; t=$(curl -s -m 60 --retry 3 "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$s&targ=self&form=text&view=full")
  st=$(echo "$t" | grep -m1 library_strategy | sed 's/.*= //'); sel=$(echo "$t" | grep -m1 library_selection | sed 's/.*= //'); mol=$(echo "$t" | grep -m1 molecule_ch1 | sed 's/.*= //' | tr -d '\r ')
  kw=$(echo "$t" | grep -i -o -E "ribo-?zero|rrna[ -]?(depl|remov)[a-z]*|ribo-?off|rnase r|poly\(?a\)?|oligo ?\(?dt\)?|truseq[^.,;]*|ribominus|total rna" | sort -u | tr '\n' ';' | cut -c1-110)
  echo -e "$g\t$s\t$st\t$sel\t$mol\t$kw"
done
