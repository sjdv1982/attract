website='http://localhost/cgi/server/ATTRACT'
rm -rf cgi-local
mkdir cgi-local
for i in `ls -d ./* | grep -v 'cgi-local' | grep -v 'serverconfig' | grep -v 'attractserver' | grep -v 'upload.py'`; do
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
cd cgi-local
ln ../upload.py
cd ..

rm -rf html-local
mkdir html-local
for i in `ls -ad html/* | \
grep -v 'full.html' | grep -v 'easy.html' | grep -v 'cryo.html' | grep -v 'upload.html' | grep -v 'peptide.html' | \
grep -v 'index.html' | grep -v 'attract.html' | grep -v 'narefine.html' | grep -v 'serverconfig' \
`; do
  cd html-local
  ln -s ../$i ${i##*/}
  cd ..
done  
cd html-local
ln -s ../html/upload-local.html upload.html
cd ..
\cp serverconfig-local.py html-local/serverconfig.py
\cp -d html/index.html html/attract.html html-local/
python generate-html-full.py $website/attractserver.py > html-local/full.html
python generate-html-easy.py $website/attractserver-easy.py > html-local/easy.html 
python generate-html-peptide.py $website/attractserver-peptide.py > html-local/peptide.html 
python generate-html-cryo.py $website/attractserver-cryo.py > html-local/cryo.html 
python generate-html-narefine.py $website/attractserver-narefine.py > html-local/narefine.html 
