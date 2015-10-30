git init
git add README.md achieve.list *.txt *.py *.sh
git commit -m "first commit"
git remote add origin https://github.com/JinfengChen/RIL_update.git
git push -u origin master


echo "update"
git add README.md achieve.list *.txt *.py *.sh
git commit -m "update"
git remote add origin https://github.com/JinfengChen/RIL_update.git
git push -u origin master
