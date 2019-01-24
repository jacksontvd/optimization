#  define a function which for a given Z and A and possibly reaction type, returns the string of the isotope, as well as the line of the parameter file which should be altered. 
def isotope(Z, A, reac_type = None):
    if int(A) == 234:
        name = 'U234nf'
        i = 1
        reaction_type = 'nf'
        reaction_word = 'induced'
        typed_name = '^{234}$U(nf)'
    elif int(A) == 236 or int(A) == 235:
        name = 'U236nf'
        i = 2
        reaction_type = 'nf'
        reaction_word = 'induced'
        typed_name = '^{236}$U(nf)'
    elif int(A) == 238 and int(Z) == 92 and reac_type == 'induced':
        name = 'U239nf'
        i = 4
        reaction_type = 'nf'
        reaction_word = 'induced'
        typed_name = '^{238}$U(nf)'
    elif int(A) == 238 and int(Z) == 92:
        name = 'U238sf'
        i = 3
        reaction_type = 'sf'
        reaction_word = 'spontaneous'
        typed_name = '^{238}$U(sf)'
    elif int(A) == 239 and int(Z) == 92:
        name = 'U239nf'
        i = 4
        reaction_type = 'nf'
        reaction_word = 'induced'
        typed_name = '^{239}$U(nf)'
    elif int(A) == 238 and int(Z) == 94:
        name = 'Pu238sf'
        i = 5
        reaction_type = 'sf'
        reaction_word = 'spontaneous'
        typed_name = '^{238}$Pu(sf)'
    elif int(A) == 239 and reac_type == 'induced':
        name = 'Pu240nf'
        i = 6
        reaction_type = 'nf'
        reaction_word = 'induced'
        typed_name = '^{240}$Pu(nf)'
    elif int(A) == 240: 
        name = 'Pu240sf'
        i = 7
        reaction_type = 'sf'
        reaction_word = 'spontaneous'
        typed_name = '^{240}$Pu(sf)'
    elif int(A) == 242 and reac_type == 'induced':
        name = 'Pu242nf'
        i = 8
        reaction_type = 'nf'
        reaction_word = 'induced'
        typed_name = '^{242}$Pu(nf)'
    elif int(A) == 242 and reac_type == 'spontaneous':
        name = 'Pu242sf'
        i = 9
        reaction_type = 'sf'
        reaction_word = 'spontaneous'
        typed_name = '^{242}$Pu(sf)'
    elif int(A) == 244:
        name =  'Cm244sf'
        i = 10
        reaction_type = 'sf'
        reaction_word = 'spontaneous'
        typed_name = '^{244}$Cm(sf)'
    elif int(A) == 252:
        name =  'Cf252sf'
        i = 11
        reaction_type = 'sf'
        reaction_word = 'spontaneous'
        typed_name = '^{252}$Cf(sf)'
    elif int(A) == 241:
        name = 'Pu242nf'
        i = 8
        reaction_type = 'nf'
        reaction_word = 'induced'
        typed_name = '^{242}$Pu(nf)'
    else:
        i = None
        name = None
        print('ERROR: Isotope either misread, or not supported by FREYA',A,A,A,A)
    return name, i, reaction_type, reaction_word, typed_name

