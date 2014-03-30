

//struct multiEpochVRcalcPSInputType
//{
//	//for both
//	vector<double> epochs;
//	bool verbose;
//	double Mass1; //in Msun
//	//for planet
//	double e_p;
//	double T_p;
//	double period_p; //in [yrs]
//	double Mass2sinI_p; //in Msun
//	double argPeri_deg_p;
//	//for companion star
//	double e_s;
//	double T_s;
//	double period_s; //in [yrs]
//	double Mass2_s; //in Msun
//	double argPeri_deg_s;
//	double inclination_deg_s;
//};

//double diff(double a, double b)
//{
//	/*
//	 This function just quickly calculates the difference between two numbers taking negatives into account.
//	 This does not use the abs() function as it is always returning 0.0 for some reason...??? don't know why.
//	 */
//
//	double difference;
//	// since the numerator is the difference, we need to take sign of each into account
//	if ( ((a>=0.0)&&(b>=0.0)) || ((a<0.0)&&(b<0.0)) )
//		// same sign so just subtract, squaring later will clear any resulting negative
//		if (a>b)
//			difference = a - b;
//		else
//			difference = b - a;
//	else if ( (a>=0.0)&&(b<=0.0) )
//		// real is negative, so subtract it
//		difference = a - b; //NOTE: tried using abs() like in Python, but it causes result to =0.0 in all cases.
//	else if ( (a<0.0)&&(b>0.0) )
//		// model is negative, so subtract it
//		difference = b - a;
//	//cout << "In diff: "<< "a ="<<a << ", b ="<<b << ", difference ="<<difference<<endl;
//	return difference;
//}
