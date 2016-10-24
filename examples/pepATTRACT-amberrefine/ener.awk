BEGIN{ifl=0}
{if($1 == "FINAL") ifl=1;
if(ifl == 1 && $1 == "BOND") {eb=$3;ea=$6;ed=$9};
if(ifl == 1 && $1 == "VDWAALS") {ev=$3;ee=$6;eg=$9};
if(ifl == 1 && $1 == "1-4") {ev14=$4;ee14=$8}}
END{printf("%13.4f%13.4f%13.4f%13.4f%13.4f%13.4f\n",
eb+ea+ed+ev+ev14+ee+ee14+eg,eb+ea+ed,ev+ev14+ee+ee14+eg,ee+ee14,eg,ee+ee14+eg)}
