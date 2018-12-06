website='http://www.attract.ph.tum.de/cgi/services/ATTRACT-devel'
python2 generate-html-full.py $website/attractserver.py > html/full.html
python2 generate-html-easy.py $website/attractserver-easy.py > html/easy.html
python2 generate-html-standard.py $website/attractserver-standard.py > html/standard.html
python2 generate-html-peptide.py $website/attractserver-peptide.py > html/peptide.html
python2 generate-html-cryo.py $website/attractserver-cryo.py > html/cryo.html
python2 generate-html-cryo-easy.py $website/attractserver-cryo-easy.py > html/cryo_easy.html 
\rm serverconfig.py
ln -s serverconfig-tum-devel.py serverconfig.py
\cp html/upload-devel.html html/upload.html
