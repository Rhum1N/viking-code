import sys




if __name__ == '__main__':
    print(sys.argv)
    a = int(sys.argv[1])
    b = int(sys.argv[2])
    #Start of the script
    #For all results file, get the actions in target.txt and draw the multiple paths.
    
    file = open("targets.txt","w")
    for i in range (a,b+1) :
        filename =  "result"+str(i)+r"/target.txt"
        f = open(filename,"r",encoding="utf-8")
        j = 0
        for ligne in f :
            file.write(ligne.rstrip() + ", ")
            j+=1
            file.write("\n")
    f.close()
    