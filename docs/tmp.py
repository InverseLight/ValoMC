from bs4 import BeautifulSoup
import codecs

with codecs.open('tmp.html','r',encoding='utf8') as file:
	file = open('tmp.html', 'r')
	str = file.read()
	file.close();

soup = BeautifulSoup(str,"lxml")
tmp=soup.find("div", {"class": "content"})
print(tmp)
