function ru=ConstructionRu(ModePropreNorme,donneeS,nbmode);

    K_elem=[1 -1;-1 1];			%matrice elementaire
	%toto.Ru=zeros(size(ModePropreNorme.Valeur,2)*(donneeS.nelem+1),donneeS.nelem);		%initialisation de la matrice Rw
   
    %initialisation de la matrice Ru
    for j=1:nbmode
        ru{j}=zeros(donneeS.nelem+1,donneeS.nelem);
    end
    
	%assemblage
    for j=1:nbmode
        for m=1:(donneeS.nelem+1)
            for k=1:donneeS.nelem
                coeff=0; 
                for l=1:donneeS.nelem
                   if l==j
                       coeff = coeff - (1/2)*ModePropreNorme.Vecteur{l}(m)*(ModePropreNorme.Vecteur{j}(k:k+1)*(((donneeS.Elem{k}.young*donneeS.Elem{k}.S/donneeS.Elem{k}.dx)*K_elem)*ModePropreNorme.Vecteur{l}(k:k+1)'));
                   else
                       coeff = coeff + ((ModePropreNorme.Valeur(l))/(ModePropreNorme.Valeur(j)-ModePropreNorme.Valeur(l)))*ModePropreNorme.Vecteur{l}(m)*(ModePropreNorme.Vecteur{j}(k:k+1)*(((donneeS.Elem{k}.young*donneeS.Elem{k}.S/donneeS.Elem{k}.dx)*K_elem)*ModePropreNorme.Vecteur{l}(k:k+1)'));
                   end
                end
                ru{j}(m,k) = ru{j}(m,k) + coeff;
            end     
        end

    end