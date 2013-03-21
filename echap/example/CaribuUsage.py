from alinea.caribu.CaribuScene import CaribuScene


#scene from a g
g=mtg_interpreter(g) #mtginterpreter calcule la geometrie
scene = plot3d(g) # on produit la scene

# source emission = lumiere intensite 1, oriente selon vecteur (0,0 -1)
source = (1,(0,0,-1))
c_scene = CaribuScene()    
idmap = c_scene.add_Shapes(scene)    
c_scene.addSources(source)
output = c_scene.runCaribu(infinity=False)
resultat = c_scene.output_by_id(output, idmap)['Einc']


# Je vien de voir que pour faire un plot de valeur 'pechees' sur g, tu pouvais jeter un oeil sur alep/test/tesprotocoladel, lignes 540 et plus

# a bientot !

#christian