template<>
EXP PosCyl toPosDeriv(const PosUVSph& p, PosDerivT<UVSph, Cyl>* deriv, PosDeriv2T<UVSph, Cyl>* deriv2)
{
	const double sn = sin(p.v), cs = cos(p.v);
	const double sh = sinh(p.u), ch = cosh(p.u);
	const double R = p.coordsys.Delta * sh * sn;
	const double z = p.coordsys.Delta * ch * cs;
	if(deriv!=NULL) {
		deriv->dRdu = p.coordsys.Delta * ch * sn;
		deriv->dRdv = p.coordsys.Delta * sh * cs;
		deriv->dzdu = p.coordsys.Delta * sh * cs;
		deriv->dzdv =-p.coordsys.Delta * ch * sn;
	}
	if(deriv2!=NULL) {
		deriv2->d2Rdu2  = R;
		deriv2->d2Rdv2  =-R;
		deriv2->d2Rdudv = z;
		deriv2->d2zdu2  = z;
		deriv2->d2zdv2  =-z;
		deriv2->d2zdudv =-R;
	}
	return PosCyl(R, z, p.phi);
}

template<>
EXP PosUVSph toPosDeriv(const PosCyl& from, const UVSph& cs,
			  PosDerivT<Cyl, UVSph>* deriv, PosDeriv2T<Cyl, UVSph>* deriv2)
{
    // lambda and nu are roots "t" of equation  R^2/(t-Delta^2) + z^2/t = 1
	double R2     = pow_2(from.R), z2 = pow_2(from.z);
	double R2_z2 = (R2+z2)/cs.Delta2;
	double shu2 = .5*(R2_z2-1+sqrt(pow_2(R2_z2)+2*(R2-z2)/cs.Delta2+1));
	double shu = sqrt(shu2);
	double chu2 = 1+shu2, chu=sqrt(chu2), ch2u=chu2+shu2;
	double cosv = from.z/(cs.Delta*chu);
	double sinv = from.R/(cs.Delta*shu), cos2v = pow_2(cosv) - pow_2(sinv);;
	double u = asinh(shu);
	double v = acos(cosv);
	if(deriv!=NULL){
		deriv->dudR = from.R * chu / (shu*(cs.Delta2 * ch2u - R2 - z2));
		deriv->dudz = from.z * shu / (chu*(cs.Delta2 * ch2u - R2 - z2));
		deriv->dvdR = from.R * cosv / (sinv*(R2+z2-cs.Delta2*cos2v));
		deriv->dvdz = from.z * sinv / (cosv*(R2+z2-cs.Delta2*cos2v));
	}
	if(deriv2!=NULL) {
		//Missing code
	}
	return PosUVSph(u, v, from.phi, cs);
}

