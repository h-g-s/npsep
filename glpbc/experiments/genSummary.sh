echo "" > tresults.txt
for file in exp*;
do
   cat $file/results.txt >> tresults.txt
done 

txt2tags -t html -o results.html tresults.txt
txt2tags -t adoc -o results.txt tresults.txt 
firefox results.html 
