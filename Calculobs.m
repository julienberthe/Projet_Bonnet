function toto=Calculobs(Vecteurp,ModePropreNorme,matriceR);

toto.Valeur = ModePropreNorme.Valeur' + matriceR.Rw*Vecteurp;

for i=1:size(ModePropreNorme.Vecteur,2)
    for j=1:size(ModePropreNorme.Vecteur{i},2)
    toto.Vecteur{i}(j) = ModePropreNorme.Vecteur{i}(j) + matriceR.Ru((i-1)*size(ModePropreNorme.Vecteur,2)+j,:)*Vecteurp;
    end
end

end