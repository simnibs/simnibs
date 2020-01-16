// simnibs automatically generated file for getdp simulation
// file created with getdp.py on 20150309_175222


WM = 1;
WM_surface = 1001;
GM = 2;
GM_surface = 1002;
CSF = 3;
CSF_surface = 1003;
Bone = 4;
Bone_surface = 1004;
Scalp = 5;
Scalp_surface = 1005;
Spongy_bone = 6;
Spongy_bone_surface = 1006;
Compact_bone = 7;
Compact_bone_surface = 1007;
Eye_balls = 8;
Eye_balls_surface = 1008;
Eye_region = 9;
Eye_region_surface = 1009;
Electrode_rubber = 100;
Electrode_rubber_surface = 1100;
Saline = 500;
Saline_surface = 1500;


fn_tensor = "vn_tensor.msh"; // tensor type: vn


Group {
	WM = Region[WM];
	WM_surface = Region[WM_surface];
	GM = Region[GM];
	GM_surface = Region[GM_surface];
	CSF = Region[CSF];
	CSF_surface = Region[CSF_surface];
	Bone = Region[Bone];
	Bone_surface = Region[Bone_surface];
	Scalp = Region[Scalp];
	Scalp_surface = Region[Scalp_surface];
	Spongy_bone = Region[Spongy_bone];
	Spongy_bone_surface = Region[Spongy_bone_surface];
	Compact_bone = Region[Compact_bone];
	Compact_bone_surface = Region[Compact_bone_surface];
	Eye_balls = Region[Eye_balls];
	Eye_balls_surface = Region[Eye_balls_surface];
	Eye_region = Region[Eye_region];
	Eye_region_surface = Region[Eye_region_surface];
	Electrode_rubber = Region[Electrode_rubber];
	Electrode_rubber_surface = Region[Electrode_rubber_surface];
	Saline = Region[Saline];
	Saline_surface = Region[Saline_surface];
	Omega = Region[{WM,GM,CSF,Bone,Scalp,Spongy_bone,Compact_bone,Eye_balls,Eye_region,Electrode_rubber,Saline}];
	Gama = Region[{WM_surface, GM_surface, CSF_surface, Bone_surface, Scalp_surface, Spongy_bone_surface, Compact_bone_surface, Eye_balls_surface, Eye_region_surface, Electrode_rubber_surface, Saline_surface}];
}


Function {
	sigma[WM] = TensorField[XYZ[]]{10} #1 ? 0.126 * #1 : Tensor[0.126, 0, 0, 0, 0.126, 0, 0, 0, 0.126]; // vn
	sigma[GM] = TensorField[XYZ[]]{10} #1 ? 0.275 * #1 : Tensor[0.275, 0, 0, 0, 0.275, 0, 0, 0, 0.275]; // vn
	sigma[CSF] = Tensor[1.654, 0, 0, 0, 1.654, 0, 0, 0, 1.654];
	sigma[Bone] = Tensor[0.01, 0, 0, 0, 0.01, 0, 0, 0, 0.01];
	sigma[Scalp] = Tensor[0.465, 0, 0, 0, 0.465, 0, 0, 0, 0.465];
	sigma[Spongy_bone] = Tensor[0.025, 0, 0, 0, 0.025, 0, 0, 0, 0.025];
	sigma[Compact_bone] = Tensor[0.008, 0, 0, 0, 0.008, 0, 0, 0, 0.008];
	sigma[Eye_balls] = Tensor[0.5, 0, 0, 0, 0.5, 0, 0, 0, 0.5];
	sigma[Eye_region] = Tensor[0.25, 0, 0, 0, 0.25, 0, 0, 0, 0.25];
	sigma[Electrode_rubber] = Tensor[0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1];
	sigma[Saline] = Tensor[1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0];
}


Function {
	dadt[] = VectorField[XYZ[]]{0};
}


Jacobian {
	{ Name Volume;  Case {{ Region Omega; Jacobian Vol;  }} }
}


Integration {{
	Name GradGrad;
	Case {{
		Type Gauss;
		Case {
			{GeoElement Tetrahedron; NumberOfPoints 1;}
		}
	}}
}}


FunctionSpace {{
	Name Hgrad_vf_Ele;
	Type Form0;
	BasisFunction {{
		Name sn;
		NameOfCoef vn;
		Function BF_Node;
		Support Region[{Omega}];
		Entity NodesOf[All];
	}}
}}


Formulation {{
	Name QS_Formulation;
	Type FemEquation;
	Quantity {{
		Name v;
		Type Local;
		NameOfSpace Hgrad_vf_Ele;
	}}
	Equation {
		Galerkin {
			[sigma[] * Dof{Grad v}, {Grad v}];
			In Omega;
			Jacobian Volume;
			Integration GradGrad;
		}
		Galerkin {
			[sigma[] * dadt[], {Grad v}];
			In Omega;
			Jacobian Volume;
			Integration GradGrad;
		}
	}
}}


Resolution {{
	Name QS;
	System {{
		Name QS_dAdt;
		NameOfFormulation QS_Formulation;
	}}
	Operation {
		fn_dadt = "dadt.msh.msh";
		GmshRead[fn_dadt,0];
		GmshRead[fn_tensor,10];
		Generate QS_dAdt;
		Solve QS_dAdt;
		SaveSolution QS_dAdt;
	}
}}


