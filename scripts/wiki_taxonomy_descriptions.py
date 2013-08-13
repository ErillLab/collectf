#from collectfapp import models
import urllib
import urllib2
from BeautifulSoup import BeautifulSoup

def get_article_paragraph(title):
    url = "http://en.wikipedia.org"
    article = urllib.quote(title)
    opener = urllib2.build_opener()
    opener.addheaders = [('User-agent', 'Mozilla/5.0')] #wikipedia needs this
    resource = opener.open(url + '/wiki/' + article)
    data = resource.read()
    resource.close()
    soup = BeautifulSoup(data)
    div =  soup.find('div', {'id':'bodyContent'})
    for x in div.contents:
        print type(x)
    
def run():
    for m in models.Taxonomy.objects.iterator():
        print m.name
        m.description = get_article_paragraph(m.name)
        m.save()

    
    

