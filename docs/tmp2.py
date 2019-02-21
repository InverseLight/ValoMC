from bs4 import BeautifulSoup
file = open('functionreference.html', 'r')
str = file.read()
file.close()
soup = BeautifulSoup(str)
for a in soup.findAll('a'):
    a['href'] = a['href'].replace("matlab:doc('", "").replace("')",'.html')

print(soup.prettify())
file = open("testfile.txt","w") 
 
file.write(soup.prettify()) 
 
file.close() 
