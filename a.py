#/usr/bin/env python
import pickle,json,random
import sys

dic = {}
dic["AAAAAAAAAAAA"] = ["1:111111"]
dic["AAAAAAAAAAAT"] = ["1:111111","1:111112"]
dic["AAAAAAAAAAAG"] = ["1:111111","1:111112","1:111113"]

def getRandomSeq():
	s = ""
	atcg = "ATCG"
	for i in range(12):
		base = atcg[random.randint(0,100)%4]
		s += base
	return s

for i in range(100000):
	s = getRandomSeq()
	dic[s] = []
	n = random.randint(1,200)
	for j in range(n):
		pos = "1:" + str(random.randint(0,9999999))
		dic[s].append(pos)		


with open("a.pkl",'wb') as f:
	pickle.dump(dic,f)

w = open("a.json",'w')
w.write(json.dumps(dic))
w.close()

