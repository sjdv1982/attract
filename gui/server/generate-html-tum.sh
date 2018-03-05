website='http://www.attract.ph.tum.de/cgi/services/ATTRACT'
python generate-html-full.py $website/attractserver.py > html/full.html
python generate-html-easy.py $website/attractserver-easy.py > html/easy.html
python generate-html-standard.py $website/attractserver-standard.py > html/standard.html
python generate-html-peptide.py $website/attractserver-peptide.py > html/peptide.html
python generate-html-cryo.py $website/attractserver-cryo.py > html/cryo.html
python generate-html-cryo-easy.py $website/attractserver-cryo-easy.py > html/cryo_easy.html
\rm serverconfig.py
ln -s serverconfig-tum.py serverconfig.py
\cp html/upload-default.html html/upload.html
