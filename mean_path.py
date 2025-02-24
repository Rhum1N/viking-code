import sys




if __name__ == '__main__':
    print(sys.argv)
    a = int(sys.argv[1])
    b = int(sys.argv[2])
    #Start of the script
    #For all results file, get the actions in action.txt and draw the multiple paths.
    paths = []
    file = open("actions.txt","w")
    for i in range (a,b+1) :
        l = []
        filename =  "result"+str(i)+r"/action_file.txt"
        f = open(filename,"r",encoding="utf-8")
        j = 0
        for ligne in f :
            file.write(ligne.rstrip() + ", ")
            l.append(int(ligne.rstrip()))
            j+=1
        paths.append(l)
        file.write("\n")
    f.close()
    print(paths)