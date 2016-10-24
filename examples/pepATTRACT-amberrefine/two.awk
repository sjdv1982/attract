BEGIN{i=0;j=0;ibeg=0;irec=0}
{if(ibeg == 1) print $0;if($1 == "ENDMDL") {i++;ibeg=0;irec=0};if ($1 == "MODEL" && $2 == va1) irec=1;if($1 == "TER" && irec == 1) ibeg = 1}
