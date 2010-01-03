function rw=ConstructionRw(ModePropreNorme,donneeS,nbmode);

    K_elem=[1 -1;-1 1];			%matrice elementaire
	%toto.Rw=zeros(size(ModePropreNorme.Valeur,2),donneeS.nelem);		%initialisation de la matrice Rw
    
    %initialisation des matrices Rw
    for i=1:nbmode
        rw{i}=zeros(donneeS.nelem);
    end
    
	%assemblage
    for j=1:nbmode
        for k=1:donneeS.nelem %ici on fait une boucle sur tous les modes propres
            %donnee.Elem{j}.young*donnee.Elem{j}.S/donnee.Elem{j}.dx
            %toto.Rw(j,k) = toto.Rw(j,k) + ModePropreNorme.Valeur(j)*(ModePropreNorme.Vecteur{j}(k:k+1)*(((donneeS.Elem{k}.young*donneeS.Elem{k}.S/donneeS.Elem{k}.dx)*K_elem)*ModePropreNorme.Vecteur{j}(k:k+1)'));
            rw{j}(k,k) = rw{j}(k,k) + ModePropreNorme.Valeur(j)*(ModePropreNorme.Vecteur{j}(k:k+1)*(((donneeS.Elem{k}.young*donneeS.Elem{k}.S/donneeS.Elem{k}.dx)*K_elem)*ModePropreNorme.Vecteur{j}(k:k+1)'));
        end
    end
end