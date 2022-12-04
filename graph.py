file1 = open('facebook.txt', 'r')
Lines = file1.readlines()

m = len(Lines);
n = 0;
for t in Lines:
    a = t.split(" ");
    n = max(n, max(int(a[1]), int(a[1])));

print(str(n)+" "+str(m));

for t in Lines:
    a = t.split(" ");
    print(str(int(a[0])+1)+" "+str(int(a[1])+1));
