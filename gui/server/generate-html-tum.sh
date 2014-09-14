website='http://www.attract.ph.tum.de/cgi/services/ATTRACT'
python generate-html-full.py $website/attractserver.py > html/full.html
python generate-html-easy.py $website/attractserver-easy.py > html/easy.html 
python generate-html-cryo.py $website/attractserver-cryo.py > html/cryo.html
python generate-html-narefine.py $website/attractserver-narefine.py > html/narefine.html 
