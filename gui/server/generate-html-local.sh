website='http://localhost/cgi/server/ATTRACT'
rm -rf cgi-local
mkdir cgi-local
for i in `ls -d ./* | grep -v 'cgi-local' | grep -v 'serverconfig' | grep -v 'attractserver'`; do
 cd cgi-local
 ln -s ../$i
 cd ..
done
cd cgi-local
ln -s ../serverconfig-local.py serverconfig.py 
cd ..
for i in `ls -d ./* | grep 'attractserver'`; do
 cd cgi-local
 ln ../$i
 cd ..
done

rm -rf html-local
mkdir html-local
for i in `ls -ad html/* | grep -v 'full.html' | grep -v 'easy.html' | grep -v 'cryo.html' | grep -v 'serverconfig' `; do
  cd html-local
  ln -s ../$i ${i##*/}
  cd ..
done  
\cp serverconfig-local.py html-local/serverconfig.py
python generate-html-full.py $website/attractserver.py > html-local/full.html
python generate-html-easy.py $website/attractserver-easy.py > html-local/easy.html 
python generate-html-cryo.py $website/attractserve-cryo.py > html-local/cryo.html 
