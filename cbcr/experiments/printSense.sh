for file in ~/Dropbox/inst/miplibClique/*.lp;
do
   probName=`basename ${file} .lp`
   strMin=`cat $file | grep -i minimize | cut -d " " -f 1`
   if [ -n "${strMin}" ]; then
      echo $probName $strMin $strMax  MINI MINI
   else
      echo $probName $strMin $strMax  MAX MAX
   fi

done

