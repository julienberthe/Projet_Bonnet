function toto=NormalisationMode(ModePropre,matriceS);

for i=1:size(ModePropre.Vecteur,2)
    
    A=ModePropre.Vecteur{i}'*(matriceS.K_ef*ModePropre.Vecteur{i});
    toto.Valeur(i) = ModePropre.Valeur(i);
    toto.Vecteur{i}(1)=0;
    for j=1:size(ModePropre.Vecteur,2)
    toto.Vecteur{i}(j+1) = ModePropre.Vecteur{i}(j)/A;
    end
end
end