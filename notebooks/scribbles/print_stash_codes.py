import iris
### this was the ap4 file
fn = "bc179a.p41930oct.pp"


cl = iris.load(fn)

my_list_of_stash = ["m01s34i027"]


list_of_stash = []


for c in cl:
##c = cl[0]
##print( c.name())
    list_of_stash.append(c.attributes["STASH"])
    


### list of stash now contains alll the stash codes in the file


for stash in my_list_of_stash:
    if stash in list_of_stash:
        ### do stuff
    else:
        ## try another stream

        
