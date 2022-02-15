#!/usr/bin/awk -f
BEGIN{
    z_min=4225;z_max=5875;dz=50;z_ref=5000;
    t_min=0; t=t_min; t_max=1; dt=1/365.0;t=0;pi=3.14159;
    gamma=0.0064;T_bias_max=10;T_bias_min=-10;
    ddf_snow0=2;ddf_snow1=8;
    ddf_ice0=4;ddf_ice1=16;
}
function T_local(z,t) {return T[t]-gamma*(z-z_ref)+T_bias;} #temperature

NR==FNR{						
    z=$1;a[z]=$2;rain[z]=0;snow[z]=0;
    snow_melt[z]=0;ice_melt[z]=0;next;
}
NR!=FNR{					
    P[t]=$2/1000.;
    T[t]=$3;
    t+=dt;

}
END{
    for(ddf_snow=ddf_snow0;ddf_snow<ddf_snow1;ddf_snow+=1){
	for(ddf_ice=ddf_ice0;ddf_ice<ddf_ice1;ddf_ice+=1){
	    #fixing steady state
	    mb=10000;T1=T_bias_max;T0=T_bias_min;count=0;scaleP=1;
	    while((mb>0.01||mb<-0.01)&&count<16){		
		T_bias=.5*(T0+T1);count++;
		total_rf=0;total_sf=0;total_sm=0;total_im=0;for(z=z_max;z>=z_min;z-=dz)snow[z]=0.;
		for(z=z_max;z>=z_min;z-=dz){
		    for(t=0;t<t_max;t+=dt){
			p=P[t];
			if(T_local(z,t)<=0){snow[z]+=p;total_sf+=(p*a[z]);}
			else {
			    total_rf+=(p*a[z]);
			    if(snow[z]>0){
				sm=.001*ddf_snow*T_local(z,t);
				snow[z]-=sm;total_sm+=(sm*a[z]);
				if(snow[z]<0){total_sm+=(snow[z]*a[z]);snow[z]=0;};
			    } else total_im+=(.001*ddf_ice*T_local(z,t)*a[z]);
			};
		    };#t
		}#z
		mb=total_sf-(total_sm+total_im);
		if(mb>0)T0=T_bias;
		else T1=T_bias;
	    }#T_bias 
	    print ddf_snow,ddf_ice,T_bias,scaleP,total_rf,total_sf,total_sm,total_im,total_sf-total_sm-total_im,total_rf+total_sm+total_im
	    if(T_bias==T_bias_max||T_bias==T_bias_min) print "WARNING1"

	    for(scaleP=.9;scaleP<=1.1;scaleP+=.05){
		# pertubed P 
		total_rf=0;total_sf=0;total_sm=0;total_im=0;for(z=z_max;z>=z_min;z-=dz)snow[z]=0.;
		for(z=z_max;z>=z_min;z-=dz){
		    for(t=0;t<t_max;t+=dt){
			p=scaleP*P[t];
			if(T_local(z,t)<=0){snow[z]+=p;total_sf+=(p*a[z]);}
			else {
			    total_rf+=(p*a[z]);
			    if(snow[z]>0){
				sm=.001*ddf_snow*T_local(z,t);
				snow[z]-=sm;total_sm+=(sm*a[z]);
				if(snow[z]<0){total_sm+=(snow[z]*a[z]);snow[z]=0;};
			    } else total_im+=(.001*ddf_ice*T_local(z,t)*a[z]);
			};
		    };#t
		};#z
		print ddf_snow,ddf_ice,T_bias,scaleP,total_rf,total_sf,total_sm,total_im,total_sf-total_sm-total_im,total_rf+total_sm+total_im >>"out1.txt"; 
	    };# scaleP 
	    print "\n">>"out1.txt";
	}#ddf_ice
 }#ddf_snow
}#END
