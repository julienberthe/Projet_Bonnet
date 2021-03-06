function ModePropre=CalculModePropre(matrice,donnee)

[VecteurPropre,ValeurPropre]=eig(matrice.M^-1*matrice.K_ef);
%  keyboard
for i=1:size(ValeurPropre,1)
	vp(i)=ValeurPropre(i,i);
end

[VP,index]=sort(vp,'ascend');

ModePropre.n=size(vp,2);
ModePropre.Matrice=VecteurPropre(:,index);
	
for i=1:size(ValeurPropre,1)
	ModePropre.Vecteur{i}	=ModePropre.Matrice(:,i);
	ModePropre.Valeur(i)	=VP(i);
end

end