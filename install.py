import sys
import  os

programDirectory = os.path.dirname(os.path.abspath(__file__))
template=""
for line in open("template_conf/FindSV.config"):
    template += line

template=template.replace("BootstrapAnn=\"\"","BootstrapAnn=\"{}/BootstrapAnn/BootstrapAnn.py\"".format(programDirectory))

print template



