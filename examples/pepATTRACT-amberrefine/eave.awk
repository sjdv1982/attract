{if ($1 == "A" && $2 == "V") iflag=1;if(iflag == 1 && $1 == "EAMBER") {print $4;iflag=0}}
