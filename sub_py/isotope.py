#  define a function which for a given Z and A and possibly reaction type, returns the string of the isotope, as well as the line of the parameter file which should be altered. 
def isotope(Z, A, reac_type = None):
    if int(A) == 234:
        name = 'U234nf'
        i = 1
        reaction_type = '(n,f)'
        reaction_word = 'induced'
    elif int(A) == 236:
        name = 'U236nf'
        i = 2
        reaction_type = '(n,f)'
        reaction_word = 'induced'
    elif int(A) == 238 and int(Z) == 92 :
        name = 'U238sf'
        i = 3
        reaction_type = 'sf'
        reaction_word = 'spontaneous'
    elif int(A) == 239:
        name = 'U239nf'
        i = 4
        reaction_type = '(n,f)'
        reaction_word = 'induced'
    elif int(A) == 238 and int(Z) == 94:
        name = 'Pu238sf'
        i = 5
        reaction_type = 'sf'
        reaction_word = 'spontaneous'
    elif int(A) == 240 and reac_type == 'induced':
        name = 'Pu240sf'
        i = 6
        reaction_type = '(n,f)'
        reaction_word = 'induced'
    elif int(A) == 240: 
        name = 'Pu240sf'
        i = 7
        reaction_type = 'sf'
        reaction_word = 'spontaneous'
    elif int(A) == 242 and reac_type == 'induced':
        name = 'Pu242nf'
        i = 8
        reaction_type = '(n,f)'
        reaction_word = 'induced'
    elif int(A) == 242 and reac_type == 'spontaneous':
        name = 'Pu242sf'
        i = 9
        reaction_type = 'sf'
        reaction_word = 'spontaneous'
    elif int(A) == 244:
        name =  'Cm244sf'
        i = 10
        reaction_type = 'sf'
        reaction_word = 'spontaneous'
    elif int(A) == 252:
        name =  'Cf252sf'
        i = 11
        reaction_type = 'sf'
        reaction_word = 'spontaneous'
    else:
        i = None
        name = None
        print('ERROR: Isotope either misread, or not supported by FREYA')
    return name, i, reaction_type, reaction_word

