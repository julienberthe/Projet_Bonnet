function toto=ConstructionRu(ModePropreNorme,matrice,donnee,matriceR1);

% A COMPLETEMENT MODIFIER -> RIEN N'EST FAIT...

    K_elem=[1 -1;-1 1];			%matrice elementaire
	toto.Ru=zeros(size(ModePropreNorme.Valeur,2)*donneeS.nelem);		%initialisation de la matrice Rw
    
	%assemblage
    for j=1:size(ModePropreNorme.Valeur,2)
        coeff=0;
        for k=1:donneeS.nelem %ici on fait une boucle sur tous les modes propres
            %donnee.Elem{j}.young*donnee.Elem{j}.S/donnee.Elem{j}.dx
            toto.Rw((j-1)*donneeS.nelem+k,(j-1)*donneeS.nelem+k) = toto.Rw((j-1)*donneeS.nelem+k,(j-1)*donneeS.nelem+k) + ModePropreNorme.Valeur(j)*(ModePropreNorme.Vecteur{j}(k:k+1)*(((donneeS.Elem{k}.young*donneeS.Elem{k}.S/donneeS.Elem{k}.dx)*K_elem)*ModePropreNorme.Vecteur{j}(k:k+1)'));
        end
    end

end