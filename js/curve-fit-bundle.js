(function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(require,module,exports){
/**
 * Created by acastillo on 8/5/15.
 */
var Matrix = require("ml-matrix");
var math = require("./algebra");

var DEBUG = false;
/** Levenberg Marquardt curve-fitting: minimize sum of weighted squared residuals
 ----------  INPUT  VARIABLES  -----------
 func   = function of n independent variables, 't', and m parameters, 'p',
 returning the simulated model: y_hat = func(t,p,c)
 p      = n-vector of initial guess of parameter values
 t      = m-vectors or matrix of independent variables (used as arg to func)
 y_dat  = m-vectors or matrix of data to be fit by func(t,p)
 weight = weighting vector for least squares fit ( weight >= 0 ) ...
 inverse of the standard measurement errors
 Default:  sqrt(d.o.f. / ( y_dat' * y_dat ))
 dp     = fractional increment of 'p' for numerical derivatives
 dp(j)>0 central differences calculated
 dp(j)<0 one sided 'backwards' differences calculated
 dp(j)=0 sets corresponding partials to zero; i.e. holds p(j) fixed
 Default:  0.001;
 p_min  = n-vector of lower bounds for parameter values
 p_max  = n-vector of upper bounds for parameter values
 c      = an optional matrix of values passed to func(t,p,c)
 opts   = vector of algorithmic parameters
 parameter    defaults    meaning
 opts(1)  =  prnt            3        >1 intermediate results; >2 plots
 opts(2)  =  MaxIter      10*Npar     maximum number of iterations
 opts(3)  =  epsilon_1       1e-3     convergence tolerance for gradient
 opts(4)  =  epsilon_2       1e-3     convergence tolerance for parameters
 opts(5)  =  epsilon_3       1e-3     convergence tolerance for Chi-square
 opts(6)  =  epsilon_4       1e-2     determines acceptance of a L-M step
 opts(7)  =  lambda_0        1e-2     initial value of L-M paramter
 opts(8)  =  lambda_UP_fac   11       factor for increasing lambda
 opts(9)  =  lambda_DN_fac    9       factor for decreasing lambda
 opts(10) =  Update_Type      1       1: Levenberg-Marquardt lambda update
 2: Quadratic update
 3: Nielsen's lambda update equations

 ----------  OUTPUT  VARIABLES  -----------
 p       = least-squares optimal estimate of the parameter values
 X2      = Chi squared criteria
 sigma_p = asymptotic standard error of the parameters
 sigma_y = asymptotic standard error of the curve-fit
 corr    = correlation matrix of the parameters
 R_sq    = R-squared cofficient of multiple determination
 cvg_hst = convergence history

 Henri Gavin, Dept. Civil & Environ. Engineering, Duke Univ. 22 Sep 2013
 modified from: http://octave.sourceforge.net/optim/function/leasqr.html
 using references by
 Press, et al., Numerical Recipes, Cambridge Univ. Press, 1992, Chapter 15.
 Sam Roweis       http://www.cs.toronto.edu/~roweis/notes/lm.pdf
 Manolis Lourakis http://www.ics.forth.gr/~lourakis/levmar/levmar.pdf
 Hans Nielson     http://www2.imm.dtu.dk/~hbn/publ/TR9905.ps
 Mathworks        optimization toolbox reference manual
 K. Madsen, H.B., Nielsen, and O. Tingleff
 http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf
 */
var LM = {

    optimize: function(func,p,t,y_dat,weight,dp,p_min,p_max,c,opts){

        var tensor_parameter = 0;			// set to 1 of parameter is a tensor

        var iteration  = 0;			// iteration counter
        //func_calls = 0;			// running count of function evaluations

        if((typeof p[0])!="object"){
            for(var i=0;i< p.length;i++){
                p[i]=[p[i]];
            }

        }
        //p = p(:); y_dat = y_dat(:); 		// make column vectors
        var i,k;
        var eps = 2^-52;
        var Npar   = p.length;//length(p); 			// number of parameters
        var Npnt   = y_dat.length;//length(y_dat);		// number of data points
        var p_old  = Matrix.zeros(Npar,1);		// previous set of parameters
        var y_old  = Matrix.zeros(Npnt,1);		// previous model, y_old = y_hat(t;p_old)
        var X2     = 1e-2/eps;			// a really big initial Chi-sq value
        var X2_old = 1e-2/eps;			// a really big initial Chi-sq value
        var J =  Matrix.zeros(Npnt,Npar);


        if (t.length != y_dat.length) {
            console.log('lm.m error: the length of t must equal the length of y_dat');

            length_t = t.length;
            length_y_dat = y_dat.length;
            var X2 = 0, corr = 0, sigma_p = 0, sigma_y = 0, R_sq = 0, cvg_hist = 0;
            if (!tensor_parameter) {
                return;
            }
        }

        weight = weight||Math.sqrt((Npnt-Npar+1)/(math.multiply(math.transpose(y_dat),y_dat)));
        dp = dp || 0.001;
        p_min   = p_min || math.multiply(Math.abs(p),-100);
        p_max   = p_max || math.multiply(Math.abs(p),100);
        c = c || 1;
        // Algorithmic Paramters
        //prnt MaxIter  eps1  eps2  epx3  eps4  lam0  lamUP lamDN UpdateType
        opts = opts ||[  3,10*Npar, 1e-3, 1e-3, 1e-3, 1e-2, 1e-2, 11, 9, 1 ];

        var prnt          = opts[0];	// >1 intermediate results; >2 plots
        var MaxIter       = opts[1];	// maximum number of iterations
        var epsilon_1     = opts[2];	// convergence tolerance for gradient
        var epsilon_2     = opts[3];	// convergence tolerance for parameter
        var epsilon_3     = opts[4];	// convergence tolerance for Chi-square
        var epsilon_4     = opts[5];	// determines acceptance of a L-M step
        var lambda_0      = opts[6];	// initial value of damping paramter, lambda
        var lambda_UP_fac = opts[7];	// factor for increasing lambda
        var lambda_DN_fac = opts[8];	// factor for decreasing lambda
        var Update_Type   = opts[9];	// 1: Levenberg-Marquardt lambda update
        // 2: Quadratic update
        // 3: Nielsen's lambda update equations

        if ( tensor_parameter && prnt == 3 ) prnt = 2;


        if(!dp.length || dp.length == 1){
            var dp_array = new Array(Npar);
            for(var i=0;i<Npar;i++)
                dp_array[i]=[dp];
            dp=dp_array;
        }

        // indices of the parameters to be fit
        var idx   = [];
        for(i=0;i<dp.length;i++){
            if(dp[i][0]!=0){
                idx.push(i);
            }
        }

        var Nfit = idx.length;			// number of parameters to fit
        var stop = false;				// termination flag

        var weight_sq = null;
        //console.log(weight);
        if ( !weight.length || weight.length < Npnt )	{
            // squared weighting vector
            //weight_sq = ( weight(1)*ones(Npnt,1) ).^2;
            //console.log("weight[0] "+typeof weight[0]);
            var tmp = math.multiply(Matrix.ones(Npnt,1),weight[0]);
            weight_sq = math.dotMultiply(tmp,tmp);
        }
        else{
            //weight_sq = (weight(:)).^2;
            weight_sq = math.dotMultiply(weight,weight);
        }


        // initialize Jacobian with finite difference calculation
        //console.log("J "+weight_sq);
        var result = this.lm_matx(func,t,p_old,y_old,1,J,p,y_dat,weight_sq,dp,c);
        var JtWJ = result.JtWJ,JtWdy=result.JtWdy,X2=result.Chi_sq,y_hat=result.y_hat,J=result.J;
        //[JtWJ,JtWdy,X2,y_hat,J] = this.lm_matx(func,t,p_old,y_old,1,J,p,y_dat,weight_sq,dp,c);
        //console.log(JtWJ);

        if ( Math.max(Math.abs(JtWdy)) < epsilon_1 ){
            console.log(' *** Your Initial Guess is Extremely Close to Optimal ***')
            console.log(' *** epsilon_1 = ', epsilon_1);
            stop = true;
        }


        switch(Update_Type){
            case 1: // Marquardt: init'l lambda
                lambda  = lambda_0;
                break;
            default:    // Quadratic and Nielsen
                lambda  = lambda_0 * Math.max(math.diag(JtWJ));
                nu=2;
        }
        //console.log(X2);
        X2_old = X2; // previous value of X2
        //console.log(MaxIter+" "+Npar);
        //var cvg_hst = Matrix.ones(MaxIter,Npar+3);		// initialize convergence history
        var h = null;
        while ( !stop && iteration <= MaxIter ) {		// --- Main Loop
            iteration = iteration + 1;
            // incremental change in parameters
            switch(Update_Type){
                case 1:					// Marquardt
                    //h = ( JtWJ + lambda * math.diag(math.diag(JtWJ)) ) \ JtWdy;
                    //h = math.multiply(math.inv(JtWdy),math.add(JtWJ,math.multiply(lambda,math.diag(math.diag(Npar)))));
                    h = math.solve(math.add(JtWJ,math.multiply(math.diag(math.diag(JtWJ)),lambda)),JtWdy);
                    break;
                default:					// Quadratic and Nielsen
                    //h = ( JtWJ + lambda * math.eye(Npar) ) \ JtWdy;

                    h = math.solve(math.add(JtWJ,math.multiply( Matrix.eye(Npar),lambda)),JtWdy);
            }

            /*for(var k=0;k< h.length;k++){
             h[k]=[h[k]];
             }*/
            //console.log("h "+h);
            //h=math.matrix(h);
            //  big = max(abs(h./p)) > 2;
            //this is a big step
            // --- Are parameters [p+h] much better than [p] ?
            var hidx = new Array(idx.length);
            for(k=0;k<idx.length;k++){
                if (typeof h === 'number') {
                   hidx[k] = h;
                } else {
                   hidx[k]=h[idx[k]];
                }
            }
            var p_try;
            if (p.columns == 1 && p.rows == 1 && hidx.length == 1) {
                p_try = new Matrix(1,1);
                p_try[0][0] = p[0][0] + hidx[0];
            }else
                p_try = math.add(p, hidx);// update the [idx] elements

            for(k=0;k<p_try.length;k++){
                p_try[k][0]=Math.min(Math.max(p_min[k][0],p_try[k][0]),p_max[k][0]);
            }
            // p_try = Math.min(Math.max(p_min,p_try),p_max);           // apply constraints

            var delta_y = math.subtract(y_dat, func(t,p_try,c));       // residual error using p_try
            //func_calls = func_calls + 1;
            //X2_try = delta_y' * ( delta_y .* weight_sq );  // Chi-squared error criteria

            var X2_try = math.multiply(math.transpose(delta_y),math.dotMultiply(delta_y,weight_sq));

            if ( Update_Type == 2 ){  			  // Quadratic
                //    One step of quadratic line update in the h direction for minimum X2
                //var alpha =  JtWdy'*h / ( (X2_try - X2)/2 + 2*JtWdy'*h ) ;
                var JtWdy_th = math.multiply(math.transpose(JtWdy),h);
                var alpha =  math.multiply(JtWdy_th,math.inv(math.add(math.multiply(math.subtract(X2_try - X2),1/2)),math.multiply(JtWdy_th,2)));//JtWdy'*h / ( (X2_try - X2)/2 + 2*JtWdy'*h ) ;

                h = math.multiply(alpha, h);
                for(var k=0;k<idx.length;k++){
                    hidx[k]=h[idx[k]];
                }

                p_try = math.add(p ,hidx);                     // update only [idx] elements
                p_try = math.min(math.max(p_min,p_try),p_max);          // apply constraints

                delta_y = math.subtract(y_dat, func(t,p_try,c));      // residual error using p_try
                // func_calls = func_calls + 1;
                //X2_try = delta_y' * ( delta_y .* weight_sq ); // Chi-squared error criteria
                X2_try = math.multiply(math.transpose(delta_y), mat.dotMultiply(delta_y, weight_sq));
            }

            //rho = (X2 - X2_try) / ( 2*h' * (lambda * h + JtWdy) ); // Nielsen
            var rho = (X2-X2_try)/math.multiply(math.multiply(math.transpose(h),2),math.add(math.multiply(lambda, h),JtWdy));
            //console.log("rho "+rho);
            if ( rho > epsilon_4 ) {		// it IS significantly better
                //console.log("Here");
                dX2 = X2 - X2_old;
                X2_old = X2;
                p_old = p;
                y_old = y_hat;
                p = p_try;			// accept p_try

                result = this.lm_matx(func, t, p_old, y_old, dX2, J, p, y_dat, weight_sq, dp, c);
                JtWJ = result.JtWJ,JtWdy=result.JtWdy,X2=result.Chi_sq,y_hat=result.y_hat,J=result.J;
                // decrease lambda ==> Gauss-Newton method

                switch (Update_Type) {
                    case 1:							// Levenberg
                        lambda = Math.max(lambda / lambda_DN_fac, 1.e-7);
                        break;
                    case 2:							// Quadratic
                        lambda = Math.max(lambda / (1 + alpha), 1.e-7);
                        break;
                    case 3:							// Nielsen
                        lambda = math.multiply(Math.max(1 / 3, 1 - (2 * rho - 1) ^ 3),lambda);
                        nu = 2;
                        break;
                }
            }
            else {					// it IS NOT better
                X2 = X2_old;			// do not accept p_try
                if (iteration%(2 * Npar)==0) {	// rank-1 update of Jacobian
                    result = this.lm_matx(func, t, p_old, y_old, -1, J, p, y_dat, weight_sq, dp, c);
                    JtWJ = result.JtWJ,JtWdy=result.JtWdy,dX2=result.Chi_sq,y_hat=result.y_hat,J=result.J;
                }

                // increase lambda  ==> gradient descent method
                switch (Update_Type) {
                    case 1:							// Levenberg
                        lambda = Math.min(lambda * lambda_UP_fac, 1.e7);
                        break;
                    case 2:							// Quadratic
                        lambda = lambda + Math.abs((X2_try - X2) / 2 / alpha);
                        break;
                    case 3:						// Nielsen
                        lambda = lambda * nu;
                        nu = 2 * nu;
                        break;
                }
            }
        }// --- End of Main Loop

        // --- convergence achieved, find covariance and confidence intervals

        // equal weights for paramter error analysis
        weight_sq = math.multiply(math.multiply(math.transpose(delta_y),delta_y), Matrix.ones(Npnt,1));

        weight_sq.apply(function(i,j){
            weight_sq[i][j] = (Npnt-Nfit+1)/weight_sq[i][j];
        });
        //console.log(weight_sq);
        result = this.lm_matx(func,t,p_old,y_old,-1,J,p,y_dat,weight_sq,dp,c);
        JtWJ = result.JtWJ,JtWdy=result.JtWdy,X2=result.Chi_sq,y_hat=result.y_hat,J=result.J;

        /*if nargout > 2				// standard error of parameters
         covar = inv(JtWJ);
         sigma_p = sqrt(diag(covar));
         end

         if nargout > 3				// standard error of the fit
         //  sigma_y = sqrt(diag(J * covar * J'));	// slower version of below
         sigma_y = zeros(Npnt,1);
         for i=1:Npnt
         sigma_y(i) = J(i,:) * covar * J(i,:)';
         end
         sigma_y = sqrt(sigma_y);
         end

         if nargout > 4				// parameter correlation matrix
         corr = covar ./ [sigma_p*sigma_p'];
         end

         if nargout > 5				// coefficient of multiple determination
         R_sq = corrcoef([y_dat y_hat]);
         R_sq = R_sq(1,2).^2;
         end

         if nargout > 6				// convergence history
         cvg_hst = cvg_hst(1:iteration,:);
         end*/

        // endfunction  # ---------------------------------------------------------- LM

        return { p:p, X2:X2, iteration:iteration };
    },

    lm_FD_J:function(func,t,p,y,dp,c) {
        // J = lm_FD_J(func,t,p,y,{dp},{c})
        //
        // partial derivatives (Jacobian) dy/dp for use with lm.m
        // computed via Finite Differences
        // Requires n or 2n function evaluations, n = number of nonzero values of dp
        // -------- INPUT VARIABLES ---------
        // func = function of independent variables, 't', and parameters, 'p',
        //        returning the simulated model: y_hat = func(t,p,c)
        // t  = m-vector of independent variables (used as arg to func)
        // p  = n-vector of current parameter values
        // y  = func(t,p,c) n-vector initialised by user before each call to lm_FD_J
        // dp = fractional increment of p for numerical derivatives
        //      dp(j)>0 central differences calculated
        //      dp(j)<0 one sided differences calculated
        //      dp(j)=0 sets corresponding partials to zero; i.e. holds p(j) fixed
        //      Default:  0.001;
        // c  = optional vector of constants passed to y_hat = func(t,p,c)
        //---------- OUTPUT VARIABLES -------
        // J  = Jacobian Matrix J(i,j)=dy(i)/dp(j)	i=1:n; j=1:m

        //   Henri Gavin, Dept. Civil & Environ. Engineering, Duke Univ. November 2005
        //   modified from: ftp://fly.cnuce.cnr.it/pub/software/octave/leasqr/
        //   Press, et al., Numerical Recipes, Cambridge Univ. Press, 1992, Chapter 15.

        var m = y.length;			// number of data points
        var n = p.length;			// number of parameters

        dp = dp || math.multiply( Matrix.ones(n, 1), 0.001);

        var ps = p.clone();//JSON.parse(JSON.stringify(p));
        //var ps = $.extend(true, [], p);
        var J = new Matrix(m,n), del =new Array(n);         // initialize Jacobian to Zero

        for (var j = 0;j < n; j++) {
            //console.log(j+" "+dp[j]+" "+p[j]+" "+ps[j]+" "+del[j]);
            del[j] = dp[j]*(1+Math.abs(p[j][0]));  // parameter perturbation
            p[j] = [ps[j][0]+del[j]];	      // perturb parameter p(j)
            //console.log(j+" "+dp[j]+" "+p[j]+" "+ps[j]+" "+del[j]);

            if (del[j] != 0){
                y1 = func(t, p, c);
                //func_calls = func_calls + 1;
                if (dp[j][0] < 0) {		// backwards difference
                    //J(:,j) = math.dotDivide(math.subtract(y1, y),del[j]);//. / del[j];
                    //console.log(del[j]);
                    //console.log(y);
                    var column = math.dotDivide(math.subtract(y1, y),del[j]);
                    for(var k=0;k< m;k++){
                        J[k][j]=column[k][0];
                    }
                    //console.log(column);
                }
                else{
                    p[j][0] = ps[j][0] - del[j];
                    //J(:,j) = (y1 - feval(func, t, p, c)). / (2. * del[j]);
                    var column = math.dotDivide(math.subtract(y1,func(t,p,c)),2*del[j]);
                    for(var k=0;k< m;k++){
                        J[k][j]=column[k][0];
                    }

                }			// central difference, additional func call
            }

            p[j] = ps[j];		// restore p(j)

        }
        //console.log("lm_FD_J: "+ JSON.stringify(J));
        return J;

    },

    // endfunction # -------------------------------------------------- LM_FD_J
    lm_Broyden_J: function(p_old,y_old,J,p,y){
        // J = lm_Broyden_J(p_old,y_old,J,p,y)
        // carry out a rank-1 update to the Jacobian matrix using Broyden's equation
        //---------- INPUT VARIABLES -------
        // p_old = previous set of parameters
        // y_old = model evaluation at previous set of parameters, y_hat(t;p_old)
        // J  = current version of the Jacobian matrix
        // p     = current  set of parameters
        // y     = model evaluation at current  set of parameters, y_hat(t;p)
        //---------- OUTPUT VARIABLES -------
        // J = rank-1 update to Jacobian Matrix J(i,j)=dy(i)/dp(j)	i=1:n; j=1:m
        //console.log(p+" X "+ p_old)
        var h  = math.subtract(p, p_old);

        //console.log("hhh "+h);
        var h_t = math.transpose(h);
        h_t.div(math.multiply(h_t,h));

        //console.log(h_t);
        //J = J + ( y - y_old - J*h )*h' / (h'*h);	// Broyden rank-1 update eq'n
        J = math.add(J, math.multiply(math.subtract(y, math.add(y_old,math.multiply(J,h))),h_t));
        return J;
        // endfunction # ---------------------------------------------- LM_Broyden_J
    },

    lm_matx : function (func,t,p_old,y_old,dX2,J,p,y_dat,weight_sq,dp,c,iteration){
        // [JtWJ,JtWdy,Chi_sq,y_hat,J] = this.lm_matx(func,t,p_old,y_old,dX2,J,p,y_dat,weight_sq,{da},{c})
        //
        // Evaluate the linearized fitting matrix, JtWJ, and vector JtWdy,
        // and calculate the Chi-squared error function, Chi_sq
        // Used by Levenberg-Marquard algorithm, lm.m
        // -------- INPUT VARIABLES ---------
        // func   = function ofpn independent variables, p, and m parameters, p,
        //         returning the simulated model: y_hat = func(t,p,c)
        // t      = m-vectors or matrix of independent variables (used as arg to func)
        // p_old  = n-vector of previous parameter values
        // y_old  = m-vector of previous model ... y_old = y_hat(t;p_old);
        // dX2    = previous change in Chi-squared criteria
        // J   = m-by-n Jacobian of model, y_hat, with respect to parameters, p
        // p      = n-vector of current  parameter values
        // y_dat  = n-vector of data to be fit by func(t,p,c)
        // weight_sq = square of the weighting vector for least squares fit ...
        //	    inverse of the standard measurement errors
        // dp     = fractional increment of 'p' for numerical derivatives
        //          dp(j)>0 central differences calculated
        //          dp(j)<0 one sided differences calculated
        //          dp(j)=0 sets corresponding partials to zero; i.e. holds p(j) fixed
        //          Default:  0.001;
        // c      = optional vector of constants passed to y_hat = func(t,p,c)
        //---------- OUTPUT VARIABLES -------
        // JtWJ	 = linearized Hessian matrix (inverse of covariance matrix)
        // JtWdy   = linearized fitting vector
        // Chi_sq = Chi-squared criteria: weighted sum of the squared residuals WSSR
        // y_hat  = model evaluated with parameters 'p'
        // J   = m-by-n Jacobian of model, y_hat, with respect to parameters, p

        //   Henri Gavin, Dept. Civil & Environ. Engineering, Duke Univ. November 2005
        //   modified from: ftp://fly.cnuce.cnr.it/pub/software/octave/leasqr/
        //   Press, et al., Numerical Recipes, Cambridge Univ. Press, 1992, Chapter 15.


        var Npnt = y_dat.length;		// number of data points
        var Npar = p.length;		// number of parameters

        dp = dp || 0.001;


        //var JtWJ = new Matrix.zeros(Npar);
        //var JtWdy  = new Matrix.zeros(Npar,1);

        var y_hat = func(t,p,c);	// evaluate model using parameters 'p'
        //func_calls = func_calls + 1;
        //console.log(J);
        if ( (iteration%(2*Npar))==0 || dX2 > 0 ) {
            //console.log("Par");
            J = this.lm_FD_J(func, t, p, y_hat, dp, c);		// finite difference
        }
        else{
            //console.log("ImPar");
            J = this.lm_Broyden_J(p_old, y_old, J, p, y_hat); // rank-1 update
        }
        //console.log(y_dat);
        //console.log(y_hat);
        var delta_y = math.subtract(y_dat, y_hat);	// residual error between model and data
        //console.log(delta_y[0][0]);
        //console.log(delta_y.rows+" "+delta_y.columns+" "+JSON.stringify(weight_sq));
        //var Chi_sq = delta_y' * ( delta_y .* weight_sq ); 	// Chi-squared error criteria
        var Chi_sq = math.multiply(math.transpose(delta_y),math.dotMultiply(delta_y,weight_sq));
        //JtWJ  = J' * ( J .* ( weight_sq * ones(1,Npar) ) );
        var Jt = math.transpose(J);

        //console.log(weight_sq);

        var JtWJ = math.multiply(Jt, math.dotMultiply(J,math.multiply(weight_sq, Matrix.ones(1,Npar))));

        //JtWdy = J' * ( weight_sq .* delta_y );
        var JtWdy = math.multiply(Jt, math.dotMultiply(weight_sq,delta_y));


        return {JtWJ:JtWJ,JtWdy:JtWdy,Chi_sq:Chi_sq,y_hat:y_hat,J:J};
        // endfunction  # ------------------------------------------------------ LM_MATX
    }



};

module.exports = LM;
},{"./algebra":2,"ml-matrix":11}],2:[function(require,module,exports){
/**
 * Created by acastillo on 8/24/15.
 */
/**
 * Non in-place function definitions, compatible with mathjs code *
 */

'use strict';

var Matrix = require('ml-matrix');

function matrix(A,B){
    return new Matrix(A,B);
}

function ones(rows, cols){
    return Matrix.ones(rows,cols);
}

function eye(rows, cols){
    return Matrix.eye(rows, cols);
}

function zeros(rows, cols){
    return Matrix.zeros(rows, cols);
}

function random(rows, cols){
    return Matrix.rand(rows,cols);
}

function transpose(A){
    if(typeof A == 'number')
        return A;
    var result = A.clone();
    return result.transpose();
}

function add(A, B){
    if(typeof A == 'number'&&typeof B === 'number')
        return A+B;
    if(typeof A == 'number')
        return this.add(B,A);

    var result = A.clone();
    if (typeof B === 'number') { // convert B to a matrix
        B2 = new Matrix(1,1);
        B2[0][0] = B;
        B = B2;
    }
    return result.add(B);

}

function subtract(A, B){
    if(typeof A == 'number'&&typeof B === 'number')
        return A-B;
    if(typeof A == 'number')
        return this.subtract(B,A);
    var result = A.clone();
    return result.sub(B);
}

function multiply(A, B){
    if(typeof A == 'number'&&typeof B === 'number')
        return A*B;
    if(typeof A == 'number')
        return this.multiply(B,A);

    var result = A.clone();

    if(typeof B === 'number')
        result.mul(B);
    else
        result = result.mmul(B);

    if(result.rows==1&&result.columns==1)
        return result[0][0];
    else
        return result;

}

function dotMultiply(A, B){
    var result = A.clone();
    return result.mul(B);
}

function dotDivide(A, B){
    var result = A.clone();
    return result.div(B);
}

function diag(A){
    var diag = null;
    var rows = A.rows, cols = A.columns, j, r;
    //It is an array
    if(typeof cols === "undefined" && (typeof A)=='object'){
        if(A[0]&&A[0].length){
            rows = A.length;
            cols = A[0].length;
            r = Math.min(rows,cols);
            diag = Matrix.zeros(cols, cols);
            for (j = 0; j < cols; j++) {
                diag[j][j]=A[j][j];
            }
        }
        else{
            cols = A.length;
            diag = Matrix.zeros(cols, cols);
            for (j = 0; j < cols; j++) {
                diag[j][j]=A[j];
            }
        }

    }
    if(rows == 1){
        diag = Matrix.zeros(cols, cols);
        for (j = 0; j < cols; j++) {
            diag[j][j]=A[0][j];
        }
    }
    else{
        if(rows>0 && cols > 0){
            r = Math.min(rows,cols);
            diag = new Array(r);
            for (j = 0; j < r; j++) {
                diag[j] = A[j][j];
            }
        } else {
            var m = new Matrix(1,1);
            m[0][0] = A;
            return m; // in case we have a number only
        }
    }
    return diag;
}

function min(A, B){
    if(typeof A==='number' && typeof B ==='number')
        return Math.min(A,B);
    var ii = A.rows, jj = A.columns;
    var result = new Matrix(ii,jj);
    for (var i = 0; i < ii; i++) {
        for (var j = 0; j < jj; j++) {
            if (A[i][j] < B[i][j]) {
                result[i][j] = A[i][j];
            }
            else{
                result[i][j] = B[i][j];
            }
        }
    }
    return result;
}

function max(A, B){
    if(typeof A==='number' && typeof B ==='number')
        return Math.max(A,B);
    var ii = A.rows, jj = A.columns;
    var result = new Matrix(ii,jj);
    for (var i = 0; i < ii; i++) {
        for (var j = 0; j < jj; j++) {
            if (A[i][j] > B[i][j]) {
                result[i][j] = A[i][j];
            }
            else{
                result[i][j] = B[i][j];
            }
        }
    }
    return result;
}

function sqrt(A){
    if(typeof A==='number' )
        return Math.sqrt(A);
    var ii = A.rows, jj = A.columns;
    var result = new Matrix(ii,jj);
    for (var i = 0; i < ii; i++) {
        for (var j = 0; j < jj; j++) {
            result[i][j] = Math.sqrt(A[i][j]);

        }
    }
    return result;
}

function abs(A){
    if(typeof A==='number' )
        return Math.abs(A);
    var ii = A.rows, jj = A.columns;
    var result = new Matrix(ii,jj);
    for (var i = 0; i < ii; i++) {
        for (var j = 0; j < jj; j++) {
            result[i][j] = Math.abs(A[i][j]);

        }
    }
    return result;
}

function exp(A){
    if(typeof A==='number' )
        return Math.sqrt(A);
    var ii = A.rows, jj = A.columns;
    var result = new Matrix(ii,jj);
    for (var i = 0; i < ii; i++) {
        for (var j = 0; j < jj; j++) {
            result[i][j] = Math.exp(A[i][j]);
        }
    }
    return result;
}

function dotPow(A, b){
    if(typeof A==='number' )
        return Math.pow(A,b);
    //console.log(A);
    var ii = A.rows, jj = A.columns;
    var result = new Matrix(ii,jj);
    for (var i = 0; i < ii; i++) {
        for (var j = 0; j < jj; j++) {
            result[i][j] = Math.pow(A[i][j],b);
        }
    }
    return result;
}

function solve(A, B){
    if (typeof A == 'object')
       return A.solve(B);
    else
       return B/A; // if its a number
}

function inv(A){
    if(typeof A ==="number")
        return 1/A;
    return A.inverse();
}

module.exports = {
    transpose:transpose,
    add:add,
    subtract:subtract,
    multiply:multiply,
    dotMultiply:dotMultiply,
    dotDivide:dotDivide,
    diag:diag,
    min:min,
    max:max,
    solve:solve,
    inv:inv,
    sqrt:sqrt,
    exp:exp,
    dotPow:dotPow,
    abs:abs,
    matrix:matrix,
    ones:ones,
    zeros:zeros,
    random:random,
    eye:eye
};

},{"ml-matrix":11}],3:[function(require,module,exports){
/* Curve fitting using Levenberg Marquard method 

   Usage: define a function that given one data point returns the values after function evaluation
     import('ml-matrix');
     var func = function(x, k) {
        var result = new Matrix(x.rows, x.columns);
        for (var i = 0; i < x.rows; i++) 
          result[i][0] = 1/(1+(x[i][0]*k[0]));
        return result; 
     }

     fit = new curve_fit();
     // fit the data with the function using -2..2 as the valid range for the parameter k
     fit.curve_fit( data, func, [-2], [2] );
     console.log('number of function evaluations: ' + fit.nfunc);
     console.log('')
     
   Use in browser (exports cf as global variable - window.cf):
     browserify curve_fit.js -o curve-fit-bundle.js 
     
     
*/

var Matrix = require("ml-matrix");
var math = require("./algebra");
var LM = require("./LM");

// var exports = module.exports = {};
var cf = cf || {};

cf.params = [];
cf.sigma_p = [];

cf.curve_fit = function (data, func, p_min, p_max) {
    lm_func = func;

    p_min = p_min || [-5];
    p_max = p_max || [5];

    var nbPoints = data.length;
    var t = math.matrix(nbPoints, 1);//[1:Npnt]';				  // independent variable
    var y_data = math.matrix(nbPoints, 1);
    var sum = 0;
    t[0][0] = data[0][0]; y_data[0][0] = data[0][1];
    var maxY = 0;
    for (var i = 0; i < nbPoints; i++) {
        t[i][0] = data[i][0];
        y_data[i][0] = data[i][1];
        if (data[i][1] > maxY)
            maxY = data[i][1];
    }
    // weighting the values does not work - would prevent elements from becoming too large
    //for (var i = 0; i < nbPoints; i++) {
    //    y_data[i][0] /= maxY
    // }
    var weight = [nbPoints / math.sqrt(y_data.dot(y_data))];
    var opts = [3, 100, 1e-3, 1e-3, 1e-3, 1e-2, 1e-2, 11, 9, 1];
    var consts = [];                         // optional vector of constants

    var p_init = math.matrix([[0.02]]);
    var p_min = math.matrix([p_min]);
    var p_max = math.matrix([p_max]);
    var p_fit = LM.optimize(lm_func, p_init, t, y_data, weight, -0.01, p_min, p_max, consts, opts);
    //cf.sigma_p = p_fit.sigma_p;
    //cf.sigma_y = p_fit.sigma_y;
    cf.X2 = p_fit.X2;
    cf.iterations = p_fit.iteration;
    cf.params = p_fit.p;
    p_fit = p_fit.p;
}

cf.getMatrix = function(rows,columns) { return new Matrix(rows, columns); };

module.exports = cf;
if (typeof window !== 'undefined') { // we are running in node
  window.cf = cf;
}
//module.exports = function() { return cf; }


},{"./LM":1,"./algebra":2,"ml-matrix":11}],4:[function(require,module,exports){
'use strict';

var Matrix = require('../matrix');

// https://github.com/lutzroeder/Mapack/blob/master/Source/CholeskyDecomposition.cs
function CholeskyDecomposition(value) {
    if (!(this instanceof CholeskyDecomposition)) {
        return new CholeskyDecomposition(value);
    }
    value = Matrix.checkMatrix(value);
    if (!value.isSymmetric())
        throw new Error('Matrix is not symmetric');

    var a = value,
        dimension = a.rows,
        l = new Matrix(dimension, dimension),
        positiveDefinite = true,
        i, j, k;

    for (j = 0; j < dimension; j++) {
        var Lrowj = l[j];
        var d = 0;
        for (k = 0; k < j; k++) {
            var Lrowk = l[k];
            var s = 0;
            for (i = 0; i < k; i++) {
                s += Lrowk[i] * Lrowj[i];
            }
            Lrowj[k] = s = (a[j][k] - s) / l[k][k];
            d = d + s * s;
        }

        d = a[j][j] - d;

        positiveDefinite &= (d > 0);
        l[j][j] = Math.sqrt(Math.max(d, 0));
        for (k = j + 1; k < dimension; k++) {
            l[j][k] = 0;
        }
    }

    if (!positiveDefinite) {
        throw new Error('Matrix is not positive definite');
    }

    this.L = l;
}

CholeskyDecomposition.prototype = {
    get lowerTriangularMatrix() {
        return this.L;
    },
    solve: function (value) {
        value = Matrix.checkMatrix(value);

        var l = this.L,
            dimension = l.rows;

        if (value.rows !== dimension) {
            throw new Error('Matrix dimensions do not match');
        }

        var count = value.columns,
            B = value.clone(),
            i, j, k;

        for (k = 0; k < dimension; k++) {
            for (j = 0; j < count; j++) {
                for (i = 0; i < k; i++) {
                    B[k][j] -= B[i][j] * l[k][i];
                }
                B[k][j] /= l[k][k];
            }
        }

        for (k = dimension - 1; k >= 0; k--) {
            for (j = 0; j < count; j++) {
                for (i = k + 1; i < dimension; i++) {
                    B[k][j] -= B[i][j] * l[i][k];
                }
                B[k][j] /= l[k][k];
            }
        }

        return B;
    }
};

module.exports = CholeskyDecomposition;

},{"../matrix":12}],5:[function(require,module,exports){
'use strict';

var Matrix = require('../matrix');
var util = require('./util');
var hypotenuse = util.hypotenuse;
var getFilled2DArray = util.getFilled2DArray;

// https://github.com/lutzroeder/Mapack/blob/master/Source/EigenvalueDecomposition.cs
function EigenvalueDecomposition(matrix) {
    if (!(this instanceof EigenvalueDecomposition)) {
        return new EigenvalueDecomposition(matrix);
    }
    matrix = Matrix.checkMatrix(matrix);
    if (!matrix.isSquare()) {
        throw new Error('Matrix is not a square matrix');
    }

    var n = matrix.columns,
        V = getFilled2DArray(n, n, 0),
        d = new Array(n),
        e = new Array(n),
        value = matrix,
        i, j;

    if (matrix.isSymmetric()) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                V[i][j] = value[i][j];
            }
        }
        tred2(n, e, d, V);
        tql2(n, e, d, V);
    }
    else {
        var H = getFilled2DArray(n, n, 0),
            ort = new Array(n);
        for (j = 0; j < n; j++) {
            for (i = 0; i < n; i++) {
                H[i][j] = value[i][j];
            }
        }
        orthes(n, H, ort, V);
        hqr2(n, e, d, V, H);
    }

    this.n = n;
    this.e = e;
    this.d = d;
    this.V = V;
}

EigenvalueDecomposition.prototype = {
    get realEigenvalues() {
        return this.d;
    },
    get imaginaryEigenvalues() {
        return this.e;
    },
    get eigenvectorMatrix() {
        if (!Matrix.isMatrix(this.V)) {
            this.V = new Matrix(this.V);
        }
        return this.V;
    },
    get diagonalMatrix() {
        var n = this.n,
            e = this.e,
            d = this.d,
            X = new Matrix(n, n),
            i, j;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                X[i][j] = 0;
            }
            X[i][i] = d[i];
            if (e[i] > 0) {
                X[i][i + 1] = e[i];
            }
            else if (e[i] < 0) {
                X[i][i - 1] = e[i];
            }
        }
        return X;
    }
};

function tred2(n, e, d, V) {

    var f, g, h, i, j, k,
        hh, scale;

    for (j = 0; j < n; j++) {
        d[j] = V[n - 1][j];
    }

    for (i = n - 1; i > 0; i--) {
        scale = 0;
        h = 0;
        for (k = 0; k < i; k++) {
            scale = scale + Math.abs(d[k]);
        }

        if (scale === 0) {
            e[i] = d[i - 1];
            for (j = 0; j < i; j++) {
                d[j] = V[i - 1][j];
                V[i][j] = 0;
                V[j][i] = 0;
            }
        } else {
            for (k = 0; k < i; k++) {
                d[k] /= scale;
                h += d[k] * d[k];
            }

            f = d[i - 1];
            g = Math.sqrt(h);
            if (f > 0) {
                g = -g;
            }

            e[i] = scale * g;
            h = h - f * g;
            d[i - 1] = f - g;
            for (j = 0; j < i; j++) {
                e[j] = 0;
            }

            for (j = 0; j < i; j++) {
                f = d[j];
                V[j][i] = f;
                g = e[j] + V[j][j] * f;
                for (k = j + 1; k <= i - 1; k++) {
                    g += V[k][j] * d[k];
                    e[k] += V[k][j] * f;
                }
                e[j] = g;
            }

            f = 0;
            for (j = 0; j < i; j++) {
                e[j] /= h;
                f += e[j] * d[j];
            }

            hh = f / (h + h);
            for (j = 0; j < i; j++) {
                e[j] -= hh * d[j];
            }

            for (j = 0; j < i; j++) {
                f = d[j];
                g = e[j];
                for (k = j; k <= i - 1; k++) {
                    V[k][j] -= (f * e[k] + g * d[k]);
                }
                d[j] = V[i - 1][j];
                V[i][j] = 0;
            }
        }
        d[i] = h;
    }

    for (i = 0; i < n - 1; i++) {
        V[n - 1][i] = V[i][i];
        V[i][i] = 1;
        h = d[i + 1];
        if (h !== 0) {
            for (k = 0; k <= i; k++) {
                d[k] = V[k][i + 1] / h;
            }

            for (j = 0; j <= i; j++) {
                g = 0;
                for (k = 0; k <= i; k++) {
                    g += V[k][i + 1] * V[k][j];
                }
                for (k = 0; k <= i; k++) {
                    V[k][j] -= g * d[k];
                }
            }
        }

        for (k = 0; k <= i; k++) {
            V[k][i + 1] = 0;
        }
    }

    for (j = 0; j < n; j++) {
        d[j] = V[n - 1][j];
        V[n - 1][j] = 0;
    }

    V[n - 1][n - 1] = 1;
    e[0] = 0;
}

function tql2(n, e, d, V) {

    var g, h, i, j, k, l, m, p, r,
        dl1, c, c2, c3, el1, s, s2,
        iter;

    for (i = 1; i < n; i++) {
        e[i - 1] = e[i];
    }

    e[n - 1] = 0;

    var f = 0,
        tst1 = 0,
        eps = Math.pow(2, -52);

    for (l = 0; l < n; l++) {
        tst1 = Math.max(tst1, Math.abs(d[l]) + Math.abs(e[l]));
        m = l;
        while (m < n) {
            if (Math.abs(e[m]) <= eps * tst1) {
                break;
            }
            m++;
        }

        if (m > l) {
            iter = 0;
            do {
                iter = iter + 1;

                g = d[l];
                p = (d[l + 1] - g) / (2 * e[l]);
                r = hypotenuse(p, 1);
                if (p < 0) {
                    r = -r;
                }

                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                dl1 = d[l + 1];
                h = g - d[l];
                for (i = l + 2; i < n; i++) {
                    d[i] -= h;
                }

                f = f + h;

                p = d[m];
                c = 1;
                c2 = c;
                c3 = c;
                el1 = e[l + 1];
                s = 0;
                s2 = 0;
                for (i = m - 1; i >= l; i--) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = hypotenuse(p, e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);

                    for (k = 0; k < n; k++) {
                        h = V[k][i + 1];
                        V[k][i + 1] = s * V[k][i] + c * h;
                        V[k][i] = c * V[k][i] - s * h;
                    }
                }

                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;

            }
            while (Math.abs(e[l]) > eps * tst1);
        }
        d[l] = d[l] + f;
        e[l] = 0;
    }

    for (i = 0; i < n - 1; i++) {
        k = i;
        p = d[i];
        for (j = i + 1; j < n; j++) {
            if (d[j] < p) {
                k = j;
                p = d[j];
            }
        }

        if (k !== i) {
            d[k] = d[i];
            d[i] = p;
            for (j = 0; j < n; j++) {
                p = V[j][i];
                V[j][i] = V[j][k];
                V[j][k] = p;
            }
        }
    }
}

function orthes(n, H, ort, V) {

    var low = 0,
        high = n - 1,
        f, g, h, i, j, m,
        scale;

    for (m = low + 1; m <= high - 1; m++) {
        scale = 0;
        for (i = m; i <= high; i++) {
            scale = scale + Math.abs(H[i][m - 1]);
        }

        if (scale !== 0) {
            h = 0;
            for (i = high; i >= m; i--) {
                ort[i] = H[i][m - 1] / scale;
                h += ort[i] * ort[i];
            }

            g = Math.sqrt(h);
            if (ort[m] > 0) {
                g = -g;
            }

            h = h - ort[m] * g;
            ort[m] = ort[m] - g;

            for (j = m; j < n; j++) {
                f = 0;
                for (i = high; i >= m; i--) {
                    f += ort[i] * H[i][j];
                }

                f = f / h;
                for (i = m; i <= high; i++) {
                    H[i][j] -= f * ort[i];
                }
            }

            for (i = 0; i <= high; i++) {
                f = 0;
                for (j = high; j >= m; j--) {
                    f += ort[j] * H[i][j];
                }

                f = f / h;
                for (j = m; j <= high; j++) {
                    H[i][j] -= f * ort[j];
                }
            }

            ort[m] = scale * ort[m];
            H[m][m - 1] = scale * g;
        }
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            V[i][j] = (i === j ? 1 : 0);
        }
    }

    for (m = high - 1; m >= low + 1; m--) {
        if (H[m][m - 1] !== 0) {
            for (i = m + 1; i <= high; i++) {
                ort[i] = H[i][m - 1];
            }

            for (j = m; j <= high; j++) {
                g = 0;
                for (i = m; i <= high; i++) {
                    g += ort[i] * V[i][j];
                }

                g = (g / ort[m]) / H[m][m - 1];
                for (i = m; i <= high; i++) {
                    V[i][j] += g * ort[i];
                }
            }
        }
    }
}

function hqr2(nn, e, d, V, H) {
    var n = nn - 1,
        low = 0,
        high = nn - 1,
        eps = Math.pow(2, -52),
        exshift = 0,
        norm = 0,
        p = 0,
        q = 0,
        r = 0,
        s = 0,
        z = 0,
        iter = 0,
        i, j, k, l, m, t, w, x, y,
        ra, sa, vr, vi,
        notlast, cdivres;

    for (i = 0; i < nn; i++) {
        if (i < low || i > high) {
            d[i] = H[i][i];
            e[i] = 0;
        }

        for (j = Math.max(i - 1, 0); j < nn; j++) {
            norm = norm + Math.abs(H[i][j]);
        }
    }

    while (n >= low) {
        l = n;
        while (l > low) {
            s = Math.abs(H[l - 1][l - 1]) + Math.abs(H[l][l]);
            if (s === 0) {
                s = norm;
            }
            if (Math.abs(H[l][l - 1]) < eps * s) {
                break;
            }
            l--;
        }

        if (l === n) {
            H[n][n] = H[n][n] + exshift;
            d[n] = H[n][n];
            e[n] = 0;
            n--;
            iter = 0;
        } else if (l === n - 1) {
            w = H[n][n - 1] * H[n - 1][n];
            p = (H[n - 1][n - 1] - H[n][n]) / 2;
            q = p * p + w;
            z = Math.sqrt(Math.abs(q));
            H[n][n] = H[n][n] + exshift;
            H[n - 1][n - 1] = H[n - 1][n - 1] + exshift;
            x = H[n][n];

            if (q >= 0) {
                z = (p >= 0) ? (p + z) : (p - z);
                d[n - 1] = x + z;
                d[n] = d[n - 1];
                if (z !== 0) {
                    d[n] = x - w / z;
                }
                e[n - 1] = 0;
                e[n] = 0;
                x = H[n][n - 1];
                s = Math.abs(x) + Math.abs(z);
                p = x / s;
                q = z / s;
                r = Math.sqrt(p * p + q * q);
                p = p / r;
                q = q / r;

                for (j = n - 1; j < nn; j++) {
                    z = H[n - 1][j];
                    H[n - 1][j] = q * z + p * H[n][j];
                    H[n][j] = q * H[n][j] - p * z;
                }

                for (i = 0; i <= n; i++) {
                    z = H[i][n - 1];
                    H[i][n - 1] = q * z + p * H[i][n];
                    H[i][n] = q * H[i][n] - p * z;
                }

                for (i = low; i <= high; i++) {
                    z = V[i][n - 1];
                    V[i][n - 1] = q * z + p * V[i][n];
                    V[i][n] = q * V[i][n] - p * z;
                }
            } else {
                d[n - 1] = x + p;
                d[n] = x + p;
                e[n - 1] = z;
                e[n] = -z;
            }

            n = n - 2;
            iter = 0;
        } else {
            x = H[n][n];
            y = 0;
            w = 0;
            if (l < n) {
                y = H[n - 1][n - 1];
                w = H[n][n - 1] * H[n - 1][n];
            }

            if (iter === 10) {
                exshift += x;
                for (i = low; i <= n; i++) {
                    H[i][i] -= x;
                }
                s = Math.abs(H[n][n - 1]) + Math.abs(H[n - 1][n - 2]);
                x = y = 0.75 * s;
                w = -0.4375 * s * s;
            }

            if (iter === 30) {
                s = (y - x) / 2;
                s = s * s + w;
                if (s > 0) {
                    s = Math.sqrt(s);
                    if (y < x) {
                        s = -s;
                    }
                    s = x - w / ((y - x) / 2 + s);
                    for (i = low; i <= n; i++) {
                        H[i][i] -= s;
                    }
                    exshift += s;
                    x = y = w = 0.964;
                }
            }

            iter = iter + 1;

            m = n - 2;
            while (m >= l) {
                z = H[m][m];
                r = x - z;
                s = y - z;
                p = (r * s - w) / H[m + 1][m] + H[m][m + 1];
                q = H[m + 1][m + 1] - z - r - s;
                r = H[m + 2][m + 1];
                s = Math.abs(p) + Math.abs(q) + Math.abs(r);
                p = p / s;
                q = q / s;
                r = r / s;
                if (m === l) {
                    break;
                }
                if (Math.abs(H[m][m - 1]) * (Math.abs(q) + Math.abs(r)) < eps * (Math.abs(p) * (Math.abs(H[m - 1][m - 1]) + Math.abs(z) + Math.abs(H[m + 1][m + 1])))) {
                    break;
                }
                m--;
            }

            for (i = m + 2; i <= n; i++) {
                H[i][i - 2] = 0;
                if (i > m + 2) {
                    H[i][i - 3] = 0;
                }
            }

            for (k = m; k <= n - 1; k++) {
                notlast = (k !== n - 1);
                if (k !== m) {
                    p = H[k][k - 1];
                    q = H[k + 1][k - 1];
                    r = (notlast ? H[k + 2][k - 1] : 0);
                    x = Math.abs(p) + Math.abs(q) + Math.abs(r);
                    if (x !== 0) {
                        p = p / x;
                        q = q / x;
                        r = r / x;
                    }
                }

                if (x === 0) {
                    break;
                }

                s = Math.sqrt(p * p + q * q + r * r);
                if (p < 0) {
                    s = -s;
                }

                if (s !== 0) {
                    if (k !== m) {
                        H[k][k - 1] = -s * x;
                    } else if (l !== m) {
                        H[k][k - 1] = -H[k][k - 1];
                    }

                    p = p + s;
                    x = p / s;
                    y = q / s;
                    z = r / s;
                    q = q / p;
                    r = r / p;

                    for (j = k; j < nn; j++) {
                        p = H[k][j] + q * H[k + 1][j];
                        if (notlast) {
                            p = p + r * H[k + 2][j];
                            H[k + 2][j] = H[k + 2][j] - p * z;
                        }

                        H[k][j] = H[k][j] - p * x;
                        H[k + 1][j] = H[k + 1][j] - p * y;
                    }

                    for (i = 0; i <= Math.min(n, k + 3); i++) {
                        p = x * H[i][k] + y * H[i][k + 1];
                        if (notlast) {
                            p = p + z * H[i][k + 2];
                            H[i][k + 2] = H[i][k + 2] - p * r;
                        }

                        H[i][k] = H[i][k] - p;
                        H[i][k + 1] = H[i][k + 1] - p * q;
                    }

                    for (i = low; i <= high; i++) {
                        p = x * V[i][k] + y * V[i][k + 1];
                        if (notlast) {
                            p = p + z * V[i][k + 2];
                            V[i][k + 2] = V[i][k + 2] - p * r;
                        }

                        V[i][k] = V[i][k] - p;
                        V[i][k + 1] = V[i][k + 1] - p * q;
                    }
                }
            }
        }
    }

    if (norm === 0) {
        return;
    }

    for (n = nn - 1; n >= 0; n--) {
        p = d[n];
        q = e[n];

        if (q === 0) {
            l = n;
            H[n][n] = 1;
            for (i = n - 1; i >= 0; i--) {
                w = H[i][i] - p;
                r = 0;
                for (j = l; j <= n; j++) {
                    r = r + H[i][j] * H[j][n];
                }

                if (e[i] < 0) {
                    z = w;
                    s = r;
                } else {
                    l = i;
                    if (e[i] === 0) {
                        H[i][n] = (w !== 0) ? (-r / w) : (-r / (eps * norm));
                    } else {
                        x = H[i][i + 1];
                        y = H[i + 1][i];
                        q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
                        t = (x * s - z * r) / q;
                        H[i][n] = t;
                        H[i + 1][n] = (Math.abs(x) > Math.abs(z)) ? ((-r - w * t) / x) : ((-s - y * t) / z);
                    }

                    t = Math.abs(H[i][n]);
                    if ((eps * t) * t > 1) {
                        for (j = i; j <= n; j++) {
                            H[j][n] = H[j][n] / t;
                        }
                    }
                }
            }
        } else if (q < 0) {
            l = n - 1;

            if (Math.abs(H[n][n - 1]) > Math.abs(H[n - 1][n])) {
                H[n - 1][n - 1] = q / H[n][n - 1];
                H[n - 1][n] = -(H[n][n] - p) / H[n][n - 1];
            } else {
                cdivres = cdiv(0, -H[n - 1][n], H[n - 1][n - 1] - p, q);
                H[n - 1][n - 1] = cdivres[0];
                H[n - 1][n] = cdivres[1];
            }

            H[n][n - 1] = 0;
            H[n][n] = 1;
            for (i = n - 2; i >= 0; i--) {
                ra = 0;
                sa = 0;
                for (j = l; j <= n; j++) {
                    ra = ra + H[i][j] * H[j][n - 1];
                    sa = sa + H[i][j] * H[j][n];
                }

                w = H[i][i] - p;

                if (e[i] < 0) {
                    z = w;
                    r = ra;
                    s = sa;
                } else {
                    l = i;
                    if (e[i] === 0) {
                        cdivres = cdiv(-ra, -sa, w, q);
                        H[i][n - 1] = cdivres[0];
                        H[i][n] = cdivres[1];
                    } else {
                        x = H[i][i + 1];
                        y = H[i + 1][i];
                        vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
                        vi = (d[i] - p) * 2 * q;
                        if (vr === 0 && vi === 0) {
                            vr = eps * norm * (Math.abs(w) + Math.abs(q) + Math.abs(x) + Math.abs(y) + Math.abs(z));
                        }
                        cdivres = cdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi);
                        H[i][n - 1] = cdivres[0];
                        H[i][n] = cdivres[1];
                        if (Math.abs(x) > (Math.abs(z) + Math.abs(q))) {
                            H[i + 1][n - 1] = (-ra - w * H[i][n - 1] + q * H[i][n]) / x;
                            H[i + 1][n] = (-sa - w * H[i][n] - q * H[i][n - 1]) / x;
                        } else {
                            cdivres = cdiv(-r - y * H[i][n - 1], -s - y * H[i][n], z, q);
                            H[i + 1][n - 1] = cdivres[0];
                            H[i + 1][n] = cdivres[1];
                        }
                    }

                    t = Math.max(Math.abs(H[i][n - 1]), Math.abs(H[i][n]));
                    if ((eps * t) * t > 1) {
                        for (j = i; j <= n; j++) {
                            H[j][n - 1] = H[j][n - 1] / t;
                            H[j][n] = H[j][n] / t;
                        }
                    }
                }
            }
        }
    }

    for (i = 0; i < nn; i++) {
        if (i < low || i > high) {
            for (j = i; j < nn; j++) {
                V[i][j] = H[i][j];
            }
        }
    }

    for (j = nn - 1; j >= low; j--) {
        for (i = low; i <= high; i++) {
            z = 0;
            for (k = low; k <= Math.min(j, high); k++) {
                z = z + V[i][k] * H[k][j];
            }
            V[i][j] = z;
        }
    }
}

function cdiv(xr, xi, yr, yi) {
    var r, d;
    if (Math.abs(yr) > Math.abs(yi)) {
        r = yi / yr;
        d = yr + r * yi;
        return [(xr + r * xi) / d, (xi - r * xr) / d];
    }
    else {
        r = yr / yi;
        d = yi + r * yr;
        return [(r * xr + xi) / d, (r * xi - xr) / d];
    }
}

module.exports = EigenvalueDecomposition;

},{"../matrix":12,"./util":9}],6:[function(require,module,exports){
'use strict';

var Matrix = require('../matrix');

// https://github.com/lutzroeder/Mapack/blob/master/Source/LuDecomposition.cs
function LuDecomposition(matrix) {
    if (!(this instanceof LuDecomposition)) {
        return new LuDecomposition(matrix);
    }
    matrix = Matrix.checkMatrix(matrix);

    var lu = matrix.clone(),
        rows = lu.rows,
        columns = lu.columns,
        pivotVector = new Array(rows),
        pivotSign = 1,
        i, j, k, p, s, t, v,
        LUrowi, LUcolj, kmax;

    for (i = 0; i < rows; i++) {
        pivotVector[i] = i;
    }

    LUcolj = new Array(rows);

    for (j = 0; j < columns; j++) {

        for (i = 0; i < rows; i++) {
            LUcolj[i] = lu[i][j];
        }

        for (i = 0; i < rows; i++) {
            LUrowi = lu[i];
            kmax = Math.min(i, j);
            s = 0;
            for (k = 0; k < kmax; k++) {
                s += LUrowi[k] * LUcolj[k];
            }
            LUrowi[j] = LUcolj[i] -= s;
        }

        p = j;
        for (i = j + 1; i < rows; i++) {
            if (Math.abs(LUcolj[i]) > Math.abs(LUcolj[p])) {
                p = i;
            }
        }

        if (p !== j) {
            for (k = 0; k < columns; k++) {
                t = lu[p][k];
                lu[p][k] = lu[j][k];
                lu[j][k] = t;
            }

            v = pivotVector[p];
            pivotVector[p] = pivotVector[j];
            pivotVector[j] = v;

            pivotSign = -pivotSign;
        }

        if (j < rows && lu[j][j] !== 0) {
            for (i = j + 1; i < rows; i++) {
                lu[i][j] /= lu[j][j];
            }
        }
    }

    this.LU = lu;
    this.pivotVector = pivotVector;
    this.pivotSign = pivotSign;
}

LuDecomposition.prototype = {
    isSingular: function () {
        var data = this.LU,
            col = data.columns;
        for (var j = 0; j < col; j++) {
            if (data[j][j] === 0) {
                return true;
            }
        }
        return false;
    },
    get determinant() {
        var data = this.LU;
        if (!data.isSquare())
            throw new Error('Matrix must be square');
        var determinant = this.pivotSign, col = data.columns;
        for (var j = 0; j < col; j++)
            determinant *= data[j][j];
        return determinant;
    },
    get lowerTriangularMatrix() {
        var data = this.LU,
            rows = data.rows,
            columns = data.columns,
            X = new Matrix(rows, columns);
        for (var i = 0; i < rows; i++) {
            for (var j = 0; j < columns; j++) {
                if (i > j) {
                    X[i][j] = data[i][j];
                } else if (i === j) {
                    X[i][j] = 1;
                } else {
                    X[i][j] = 0;
                }
            }
        }
        return X;
    },
    get upperTriangularMatrix() {
        var data = this.LU,
            rows = data.rows,
            columns = data.columns,
            X = new Matrix(rows, columns);
        for (var i = 0; i < rows; i++) {
            for (var j = 0; j < columns; j++) {
                if (i <= j) {
                    X[i][j] = data[i][j];
                } else {
                    X[i][j] = 0;
                }
            }
        }
        return X;
    },
    get pivotPermutationVector() {
        return this.pivotVector.slice();
    },
    solve: function (value) {
        value = Matrix.checkMatrix(value);

        var lu = this.LU,
            rows = lu.rows;

        if (rows !== value.rows)
            throw new Error('Invalid matrix dimensions');
        if (this.isSingular())
            throw new Error('LU matrix is singular');

        var count = value.columns,
            X = value.subMatrixRow(this.pivotVector, 0, count - 1),
            columns = lu.columns,
            i, j, k;

        for (k = 0; k < columns; k++) {
            for (i = k + 1; i < columns; i++) {
                for (j = 0; j < count; j++) {
                    X[i][j] -= X[k][j] * lu[i][k];
                }
            }
        }
        for (k = columns - 1; k >= 0; k--) {
            for (j = 0; j < count; j++) {
                X[k][j] /= lu[k][k];
            }
            for (i = 0; i < k; i++) {
                for (j = 0; j < count; j++) {
                    X[i][j] -= X[k][j] * lu[i][k];
                }
            }
        }
        return X;
    }
};

module.exports = LuDecomposition;

},{"../matrix":12}],7:[function(require,module,exports){
'use strict';

var Matrix = require('../matrix');
var hypotenuse = require('./util').hypotenuse;

//https://github.com/lutzroeder/Mapack/blob/master/Source/QrDecomposition.cs
function QrDecomposition(value) {
    if (!(this instanceof QrDecomposition)) {
        return new QrDecomposition(value);
    }
    value = Matrix.checkMatrix(value);

    var qr = value.clone(),
        m = value.rows,
        n = value.columns,
        rdiag = new Array(n),
        i, j, k, s;

    for (k = 0; k < n; k++) {
        var nrm = 0;
        for (i = k; i < m; i++) {
            nrm = hypotenuse(nrm, qr[i][k]);
        }
        if (nrm !== 0) {
            if (qr[k][k] < 0) {
                nrm = -nrm;
            }
            for (i = k; i < m; i++) {
                qr[i][k] /= nrm;
            }
            qr[k][k] += 1;
            for (j = k + 1; j < n; j++) {
                s = 0;
                for (i = k; i < m; i++) {
                    s += qr[i][k] * qr[i][j];
                }
                s = -s / qr[k][k];
                for (i = k; i < m; i++) {
                    qr[i][j] += s * qr[i][k];
                }
            }
        }
        rdiag[k] = -nrm;
    }

    this.QR = qr;
    this.Rdiag = rdiag;
}

QrDecomposition.prototype = {
    solve: function (value) {
        value = Matrix.checkMatrix(value);

        var qr = this.QR,
            m = qr.rows;

        if (value.rows !== m)
            throw new Error('Matrix row dimensions must agree');
        if (!this.isFullRank())
            throw new Error('Matrix is rank deficient');

        var count = value.columns,
            X = value.clone(),
            n = qr.columns,
            i, j, k, s;

        for (k = 0; k < n; k++) {
            for (j = 0; j < count; j++) {
                s = 0;
                for (i = k; i < m; i++) {
                    s += qr[i][k] * X[i][j];
                }
                s = -s / qr[k][k];
                for (i = k; i < m; i++) {
                    X[i][j] += s * qr[i][k];
                }
            }
        }
        for (k = n - 1; k >= 0; k--) {
            for (j = 0; j < count; j++) {
                X[k][j] /= this.Rdiag[k];
            }
            for (i = 0; i < k; i++) {
                for (j = 0; j < count; j++) {
                    X[i][j] -= X[k][j] * qr[i][k];
                }
            }
        }

        return X.subMatrix(0, n - 1, 0, count - 1);
    },
    isFullRank: function () {
        var columns = this.QR.columns;
        for (var i = 0; i < columns; i++) {
            if (this.Rdiag[i] === 0) {
                return false;
            }
        }
        return true;
    },
    get upperTriangularMatrix() {
        var qr = this.QR,
            n = qr.columns,
            X = new Matrix(n, n),
            i, j;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (i < j) {
                    X[i][j] = qr[i][j];
                } else if (i === j) {
                    X[i][j] = this.Rdiag[i];
                } else {
                    X[i][j] = 0;
                }
            }
        }
        return X;
    },
    get orthogonalMatrix() {
        var qr = this.QR,
            rows = qr.rows,
            columns = qr.columns,
            X = new Matrix(rows, columns),
            i, j, k, s;

        for (k = columns - 1; k >= 0; k--) {
            for (i = 0; i < rows; i++) {
                X[i][k] = 0;
            }
            X[k][k] = 1;
            for (j = k; j < columns; j++) {
                if (qr[k][k] !== 0) {
                    s = 0;
                    for (i = k; i < rows; i++) {
                        s += qr[i][k] * X[i][j];
                    }

                    s = -s / qr[k][k];

                    for (i = k; i < rows; i++) {
                        X[i][j] += s * qr[i][k];
                    }
                }
            }
        }
        return X;
    }
};

module.exports = QrDecomposition;

},{"../matrix":12,"./util":9}],8:[function(require,module,exports){
'use strict';

var Matrix = require('../matrix');
var util = require('./util');
var hypotenuse = util.hypotenuse;
var getFilled2DArray = util.getFilled2DArray;

// https://github.com/lutzroeder/Mapack/blob/master/Source/SingularValueDecomposition.cs
function SingularValueDecomposition(value, options) {
    if (!(this instanceof SingularValueDecomposition)) {
        return new SingularValueDecomposition(value, options);
    }
    value = Matrix.checkMatrix(value);

    options = options || {};

    var m = value.rows,
        n = value.columns,
        nu = Math.min(m, n);

    var wantu = true, wantv = true;
    if (options.computeLeftSingularVectors === false)
        wantu = false;
    if (options.computeRightSingularVectors === false)
        wantv = false;
    var autoTranspose = options.autoTranspose === true;

    var swapped = false;
    var a;
    if (m < n) {
        if (!autoTranspose) {
            a = value.clone();
            console.warn('Computing SVD on a matrix with more columns than rows. Consider enabling autoTranspose');
        } else {
            a = value.transpose();
            m = a.rows;
            n = a.columns;
            swapped = true;
            var aux = wantu;
            wantu = wantv;
            wantv = aux;
        }
    } else {
        a = value.clone();
    }

    var s = new Array(Math.min(m + 1, n)),
        U = getFilled2DArray(m, nu, 0),
        V = getFilled2DArray(n, n, 0),
        e = new Array(n),
        work = new Array(m);

    var nct = Math.min(m - 1, n);
    var nrt = Math.max(0, Math.min(n - 2, m));

    var i, j, k, p, t, ks, f, cs, sn, max, kase,
        scale, sp, spm1, epm1, sk, ek, b, c, shift, g;

    for (k = 0, max = Math.max(nct, nrt); k < max; k++) {
        if (k < nct) {
            s[k] = 0;
            for (i = k; i < m; i++) {
                s[k] = hypotenuse(s[k], a[i][k]);
            }
            if (s[k] !== 0) {
                if (a[k][k] < 0) {
                    s[k] = -s[k];
                }
                for (i = k; i < m; i++) {
                    a[i][k] /= s[k];
                }
                a[k][k] += 1;
            }
            s[k] = -s[k];
        }

        for (j = k + 1; j < n; j++) {
            if ((k < nct) && (s[k] !== 0)) {
                t = 0;
                for (i = k; i < m; i++) {
                    t += a[i][k] * a[i][j];
                }
                t = -t / a[k][k];
                for (i = k; i < m; i++) {
                    a[i][j] += t * a[i][k];
                }
            }
            e[j] = a[k][j];
        }

        if (wantu && (k < nct)) {
            for (i = k; i < m; i++) {
                U[i][k] = a[i][k];
            }
        }

        if (k < nrt) {
            e[k] = 0;
            for (i = k + 1; i < n; i++) {
                e[k] = hypotenuse(e[k], e[i]);
            }
            if (e[k] !== 0) {
                if (e[k + 1] < 0)
                    e[k] = -e[k];
                for (i = k + 1; i < n; i++) {
                    e[i] /= e[k];
                }
                e[k + 1] += 1;
            }
            e[k] = -e[k];
            if ((k + 1 < m) && (e[k] !== 0)) {
                for (i = k + 1; i < m; i++) {
                    work[i] = 0;
                }
                for (j = k + 1; j < n; j++) {
                    for (i = k + 1; i < m; i++) {
                        work[i] += e[j] * a[i][j];
                    }
                }
                for (j = k + 1; j < n; j++) {
                    t = -e[j] / e[k + 1];
                    for (i = k + 1; i < m; i++) {
                        a[i][j] += t * work[i];
                    }
                }
            }
            if (wantv) {
                for (i = k + 1; i < n; i++) {
                    V[i][k] = e[i];
                }
            }
        }
    }

    p = Math.min(n, m + 1);
    if (nct < n) {
        s[nct] = a[nct][nct];
    }
    if (m < p) {
        s[p - 1] = 0;
    }
    if (nrt + 1 < p) {
        e[nrt] = a[nrt][p - 1];
    }
    e[p - 1] = 0;

    if (wantu) {
        for (j = nct; j < nu; j++) {
            for (i = 0; i < m; i++) {
                U[i][j] = 0;
            }
            U[j][j] = 1;
        }
        for (k = nct - 1; k >= 0; k--) {
            if (s[k] !== 0) {
                for (j = k + 1; j < nu; j++) {
                    t = 0;
                    for (i = k; i < m; i++) {
                        t += U[i][k] * U[i][j];
                    }
                    t = -t / U[k][k];
                    for (i = k; i < m; i++) {
                        U[i][j] += t * U[i][k];
                    }
                }
                for (i = k; i < m; i++) {
                    U[i][k] = -U[i][k];
                }
                U[k][k] = 1 + U[k][k];
                for (i = 0; i < k - 1; i++) {
                    U[i][k] = 0;
                }
            } else {
                for (i = 0; i < m; i++) {
                    U[i][k] = 0;
                }
                U[k][k] = 1;
            }
        }
    }

    if (wantv) {
        for (k = n - 1; k >= 0; k--) {
            if ((k < nrt) && (e[k] !== 0)) {
                for (j = k + 1; j < n; j++) {
                    t = 0;
                    for (i = k + 1; i < n; i++) {
                        t += V[i][k] * V[i][j];
                    }
                    t = -t / V[k + 1][k];
                    for (i = k + 1; i < n; i++) {
                        V[i][j] += t * V[i][k];
                    }
                }
            }
            for (i = 0; i < n; i++) {
                V[i][k] = 0;
            }
            V[k][k] = 1;
        }
    }

    var pp = p - 1,
        iter = 0,
        eps = Math.pow(2, -52);
    while (p > 0) {
        for (k = p - 2; k >= -1; k--) {
            if (k === -1) {
                break;
            }
            if (Math.abs(e[k]) <= eps * (Math.abs(s[k]) + Math.abs(s[k + 1]))) {
                e[k] = 0;
                break;
            }
        }
        if (k === p - 2) {
            kase = 4;
        } else {
            for (ks = p - 1; ks >= k; ks--) {
                if (ks === k) {
                    break;
                }
                t = (ks !== p ? Math.abs(e[ks]) : 0) + (ks !== k + 1 ? Math.abs(e[ks - 1]) : 0);
                if (Math.abs(s[ks]) <= eps * t) {
                    s[ks] = 0;
                    break;
                }
            }
            if (ks === k) {
                kase = 3;
            } else if (ks === p - 1) {
                kase = 1;
            } else {
                kase = 2;
                k = ks;
            }
        }

        k++;

        switch (kase) {
            case 1: {
                f = e[p - 2];
                e[p - 2] = 0;
                for (j = p - 2; j >= k; j--) {
                    t = hypotenuse(s[j], f);
                    cs = s[j] / t;
                    sn = f / t;
                    s[j] = t;
                    if (j !== k) {
                        f = -sn * e[j - 1];
                        e[j - 1] = cs * e[j - 1];
                    }
                    if (wantv) {
                        for (i = 0; i < n; i++) {
                            t = cs * V[i][j] + sn * V[i][p - 1];
                            V[i][p - 1] = -sn * V[i][j] + cs * V[i][p - 1];
                            V[i][j] = t;
                        }
                    }
                }
                break;
            }
            case 2 : {
                f = e[k - 1];
                e[k - 1] = 0;
                for (j = k; j < p; j++) {
                    t = hypotenuse(s[j], f);
                    cs = s[j] / t;
                    sn = f / t;
                    s[j] = t;
                    f = -sn * e[j];
                    e[j] = cs * e[j];
                    if (wantu) {
                        for (i = 0; i < m; i++) {
                            t = cs * U[i][j] + sn * U[i][k - 1];
                            U[i][k - 1] = -sn * U[i][j] + cs * U[i][k - 1];
                            U[i][j] = t;
                        }
                    }
                }
                break;
            }
            case 3 : {
                scale = Math.max(Math.max(Math.max(Math.max(Math.abs(s[p - 1]), Math.abs(s[p - 2])), Math.abs(e[p - 2])), Math.abs(s[k])), Math.abs(e[k]));
                sp = s[p - 1] / scale;
                spm1 = s[p - 2] / scale;
                epm1 = e[p - 2] / scale;
                sk = s[k] / scale;
                ek = e[k] / scale;
                b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2;
                c = (sp * epm1) * (sp * epm1);
                shift = 0;
                if ((b !== 0) || (c !== 0)) {
                    shift = Math.sqrt(b * b + c);
                    if (b < 0) {
                        shift = -shift;
                    }
                    shift = c / (b + shift);
                }
                f = (sk + sp) * (sk - sp) + shift;
                g = sk * ek;
                for (j = k; j < p - 1; j++) {
                    t = hypotenuse(f, g);
                    cs = f / t;
                    sn = g / t;
                    if (j !== k) {
                        e[j - 1] = t;
                    }
                    f = cs * s[j] + sn * e[j];
                    e[j] = cs * e[j] - sn * s[j];
                    g = sn * s[j + 1];
                    s[j + 1] = cs * s[j + 1];
                    if (wantv) {
                        for (i = 0; i < n; i++) {
                            t = cs * V[i][j] + sn * V[i][j + 1];
                            V[i][j + 1] = -sn * V[i][j] + cs * V[i][j + 1];
                            V[i][j] = t;
                        }
                    }
                    t = hypotenuse(f, g);
                    cs = f / t;
                    sn = g / t;
                    s[j] = t;
                    f = cs * e[j] + sn * s[j + 1];
                    s[j + 1] = -sn * e[j] + cs * s[j + 1];
                    g = sn * e[j + 1];
                    e[j + 1] = cs * e[j + 1];
                    if (wantu && (j < m - 1)) {
                        for (i = 0; i < m; i++) {
                            t = cs * U[i][j] + sn * U[i][j + 1];
                            U[i][j + 1] = -sn * U[i][j] + cs * U[i][j + 1];
                            U[i][j] = t;
                        }
                    }
                }
                e[p - 2] = f;
                iter = iter + 1;
                break;
            }
            case 4: {
                if (s[k] <= 0) {
                    s[k] = (s[k] < 0 ? -s[k] : 0);
                    if (wantv) {
                        for (i = 0; i <= pp; i++) {
                            V[i][k] = -V[i][k];
                        }
                    }
                }
                while (k < pp) {
                    if (s[k] >= s[k + 1]) {
                        break;
                    }
                    t = s[k];
                    s[k] = s[k + 1];
                    s[k + 1] = t;
                    if (wantv && (k < n - 1)) {
                        for (i = 0; i < n; i++) {
                            t = V[i][k + 1];
                            V[i][k + 1] = V[i][k];
                            V[i][k] = t;
                        }
                    }
                    if (wantu && (k < m - 1)) {
                        for (i = 0; i < m; i++) {
                            t = U[i][k + 1];
                            U[i][k + 1] = U[i][k];
                            U[i][k] = t;
                        }
                    }
                    k++;
                }
                iter = 0;
                p--;
                break;
            }
        }
    }

    if (swapped) {
        var tmp = V;
        V = U;
        U = tmp;
    }

    this.m = m;
    this.n = n;
    this.s = s;
    this.U = U;
    this.V = V;
}

SingularValueDecomposition.prototype = {
    get condition() {
        return this.s[0] / this.s[Math.min(this.m, this.n) - 1];
    },
    get norm2() {
        return this.s[0];
    },
    get rank() {
        var eps = Math.pow(2, -52),
            tol = Math.max(this.m, this.n) * this.s[0] * eps,
            r = 0,
            s = this.s;
        for (var i = 0, ii = s.length; i < ii; i++) {
            if (s[i] > tol) {
                r++;
            }
        }
        return r;
    },
    get diagonal() {
        return this.s;
    },
    // https://github.com/accord-net/framework/blob/development/Sources/Accord.Math/Decompositions/SingularValueDecomposition.cs
    get threshold() {
        return (Math.pow(2, -52) / 2) * Math.max(this.m, this.n) * this.s[0];
    },
    get leftSingularVectors() {
        if (!Matrix.isMatrix(this.U)) {
            this.U = new Matrix(this.U);
        }
        return this.U;
    },
    get rightSingularVectors() {
        if (!Matrix.isMatrix(this.V)) {
            this.V = new Matrix(this.V);
        }
        return this.V;
    },
    get diagonalMatrix() {
        return Matrix.diag(this.s);
    },
    solve: function (value) {

        var Y = value,
            e = this.threshold,
            scols = this.s.length,
            Ls = Matrix.zeros(scols, scols),
            i;

        for (i = 0; i < scols; i++) {
            if (Math.abs(this.s[i]) <= e) {
                Ls[i][i] = 0;
            } else {
                Ls[i][i] = 1 / this.s[i];
            }
        }

        var U = this.U;
        var V = this.rightSingularVectors;

        var VL = V.mmul(Ls),
            vrows = V.rows,
            urows = U.length,
            VLU = Matrix.zeros(vrows, urows),
            j, k, sum;

        for (i = 0; i < vrows; i++) {
            for (j = 0; j < urows; j++) {
                sum = 0;
                for (k = 0; k < scols; k++) {
                    sum += VL[i][k] * U[j][k];
                }
                VLU[i][j] = sum;
            }
        }

        return VLU.mmul(Y);
    },
    solveForDiagonal: function (value) {
        return this.solve(Matrix.diag(value));
    },
    inverse: function () {
        var V = this.V;
        var e = this.threshold,
            vrows = V.length,
            vcols = V[0].length,
            X = new Matrix(vrows, this.s.length),
            i, j;

        for (i = 0; i < vrows; i++) {
            for (j = 0; j < vcols; j++) {
                if (Math.abs(this.s[j]) > e) {
                    X[i][j] = V[i][j] / this.s[j];
                } else {
                    X[i][j] = 0;
                }
            }
        }

        var U = this.U;

        var urows = U.length,
            ucols = U[0].length,
            Y = new Matrix(vrows, urows),
            k, sum;

        for (i = 0; i < vrows; i++) {
            for (j = 0; j < urows; j++) {
                sum = 0;
                for (k = 0; k < ucols; k++) {
                    sum += X[i][k] * U[j][k];
                }
                Y[i][j] = sum;
            }
        }

        return Y;
    }
};

module.exports = SingularValueDecomposition;

},{"../matrix":12,"./util":9}],9:[function(require,module,exports){
'use strict';

exports.hypotenuse = function hypotenuse(a, b) {
    if (Math.abs(a) > Math.abs(b)) {
        var r = b / a;
        return Math.abs(a) * Math.sqrt(1 + r * r);
    }
    if (b !== 0) {
        var r = a / b;
        return Math.abs(b) * Math.sqrt(1 + r * r);
    }
    return 0;
};

// For use in the decomposition algorithms. With big matrices, access time is
// too long on elements from array subclass
// todo check when it is fixed in v8
// http://jsperf.com/access-and-write-array-subclass
exports.getEmpty2DArray = function (rows, columns) {
    var array = new Array(rows);
    for (var i = 0; i < rows; i++) {
        array[i] = new Array(columns);
    }
    return array;
};

exports.getFilled2DArray = function (rows, columns, value) {
    var array = new Array(rows);
    for (var i = 0; i < rows; i++) {
        array[i] = new Array(columns);
        for (var j = 0; j < columns; j++) {
            array[i][j] = value;
        }
    }
    return array;
};

},{}],10:[function(require,module,exports){
'use strict';

var Matrix = require('./matrix');

var SingularValueDecomposition = require('./dc/svd');
var EigenvalueDecomposition = require('./dc/evd');
var LuDecomposition = require('./dc/lu');
var QrDecomposition = require('./dc/qr');
var CholeskyDecomposition = require('./dc/cholesky');

function inverse(matrix) {
    matrix = Matrix.checkMatrix(matrix);
    return solve(matrix, Matrix.eye(matrix.rows));
}

Matrix.inverse = Matrix.inv = inverse;
Matrix.prototype.inverse = Matrix.prototype.inv = function () {
    return inverse(this);
};

function solve(leftHandSide, rightHandSide) {
    leftHandSide = Matrix.checkMatrix(leftHandSide);
    rightHandSide = Matrix.checkMatrix(rightHandSide);
    return leftHandSide.isSquare() ? new LuDecomposition(leftHandSide).solve(rightHandSide) : new QrDecomposition(leftHandSide).solve(rightHandSide);
}

Matrix.solve = solve;
Matrix.prototype.solve = function (other) {
    return solve(this, other);
};

module.exports = {
    SingularValueDecomposition: SingularValueDecomposition,
    SVD: SingularValueDecomposition,
    EigenvalueDecomposition: EigenvalueDecomposition,
    EVD: EigenvalueDecomposition,
    LuDecomposition: LuDecomposition,
    LU: LuDecomposition,
    QrDecomposition: QrDecomposition,
    QR: QrDecomposition,
    CholeskyDecomposition: CholeskyDecomposition,
    CHO: CholeskyDecomposition,
    inverse: inverse,
    solve: solve
};

},{"./dc/cholesky":4,"./dc/evd":5,"./dc/lu":6,"./dc/qr":7,"./dc/svd":8,"./matrix":12}],11:[function(require,module,exports){
'use strict';

module.exports = require('./matrix');
module.exports.Decompositions = module.exports.DC = require('./decompositions');

},{"./decompositions":10,"./matrix":12}],12:[function(require,module,exports){
'use strict';

/**
 * Real matrix
 */
class Matrix extends Array {
    /**
     * @constructor
     * @param {number|Array|Matrix} nRows - Number of rows of the new matrix,
     * 2D array containing the data or Matrix instance to clone
     * @param {number} [nColumns] - Number of columns of the new matrix
     */
    constructor(nRows, nColumns) {
        if (Matrix.isMatrix(nRows)) {
            return nRows.clone();
        } else if (Number.isInteger(nRows) && nRows > 0) { // Create an empty matrix
            super(nRows);
            if (Number.isInteger(nColumns) && nColumns > 0) {
                for (var i = 0; i < nRows; i++) {
                    this[i] = new Array(nColumns);
                }
            } else {
                throw new TypeError('nColumns must be a positive integer');
            }
        } else if (Array.isArray(nRows)) { // Copy the values from the 2D array
            var matrix = nRows;
            nRows = matrix.length;
            nColumns = matrix[0].length;
            if (typeof nColumns !== 'number' || nColumns === 0) {
                throw new TypeError('Data must be a 2D array with at least one element');
            }
            super(nRows);
            for (var i = 0; i < nRows; i++) {
                if (matrix[i].length !== nColumns) {
                    throw new RangeError('Inconsistent array dimensions');
                }
                this[i] = [].concat(matrix[i]);
            }
        } else {
            throw new TypeError('First argument must be a positive number or an array');
        }
        this.rows = nRows;
        this.columns = nColumns;
    }

    /**
     * Constructs a Matrix with the chosen dimensions from a 1D array
     * @param {number} newRows - Number of rows
     * @param {number} newColumns - Number of columns
     * @param {Array} newData - A 1D array containing data for the matrix
     * @returns {Matrix} - The new matrix
     */
    static from1DArray(newRows, newColumns, newData) {
        var length = newRows * newColumns;
        if (length !== newData.length) {
            throw new RangeError('Data length does not match given dimensions');
        }
        var newMatrix = new Matrix(newRows, newColumns);
        for (var row = 0; row < newRows; row++) {
            for (var column = 0; column < newColumns; column++) {
                newMatrix[row][column] = newData[row * newColumns + column];
            }
        }
        return newMatrix;
    }

    /**
     * Creates a row vector, a matrix with only one row.
     * @param {Array} newData - A 1D array containing data for the vector
     * @returns {Matrix} - The new matrix
     */
    static rowVector(newData) {
        var vector = new Matrix(1, newData.length);
        for (var i = 0; i < newData.length; i++) {
            vector[0][i] = newData[i];
        }
        return vector;
    }

    /**
     * Creates a column vector, a matrix with only one column.
     * @param {Array} newData - A 1D array containing data for the vector
     * @returns {Matrix} - The new matrix
     */
    static columnVector(newData) {
        var vector = new Matrix(newData.length, 1);
        for (var i = 0; i < newData.length; i++) {
            vector[i][0] = newData[i];
        }
        return vector;
    }

    /**
     * Creates an empty matrix with the given dimensions. Values will be undefined. Same as using new Matrix(rows, columns).
     * @param {number} rows - Number of rows
     * @param {number} columns - Number of columns
     * @returns {Matrix} - The new matrix
     */
    static empty(rows, columns) {
        return new Matrix(rows, columns);
    }

    /**
     * Creates a matrix with the given dimensions. Values will be set to zero.
     * @param {number} rows - Number of rows
     * @param {number} columns - Number of columns
     * @returns {Matrix} - The new matrix
     */
    static zeros(rows, columns) {
        return Matrix.empty(rows, columns).fill(0);
    }

    /**
     * Creates a matrix with the given dimensions. Values will be set to one.
     * @param {number} rows - Number of rows
     * @param {number} columns - Number of columns
     * @returns {Matrix} - The new matrix
     */
    static ones(rows, columns) {
        return Matrix.empty(rows, columns).fill(1);
    }

    /**
     * Creates a matrix with the given dimensions. Values will be randomly set.
     * @param {number} rows - Number of rows
     * @param {number} columns - Number of columns
     * @param {function} [rng] - Random number generator (default: Math.random)
     * @returns {Matrix} The new matrix
     */
    static rand(rows, columns, rng) {
        if (rng === undefined) rng = Math.random;
        var matrix = Matrix.empty(rows, columns);
        for (var i = 0; i < rows; i++) {
            for (var j = 0; j < columns; j++) {
                matrix[i][j] = rng();
            }
        }
        return matrix;
    }

    /**
     * Creates an identity matrix with the given dimension. Values of the diagonal will be 1 and others will be 0.
     * @param {number} rows - Number of rows
     * @param {number} [columns] - Number of columns (Default: rows)
     * @returns {Matrix} - The new identity matrix
     */
    static eye(rows, columns) {
        if (columns === undefined) columns = rows;
        var min = Math.min(rows, columns);
        var matrix = Matrix.zeros(rows, columns);
        for (var i = 0; i < min; i++) {
            matrix[i][i] = 1;
        }
        return matrix;
    }

    /**
     * Creates a diagonal matrix based on the given array.
     * @param {Array} data - Array containing the data for the diagonal
     * @param {number} [rows] - Number of rows (Default: data.length)
     * @param {number} [columns] - Number of columns (Default: rows)
     * @returns {Matrix} - The new diagonal matrix
     */
    static diag(data, rows, columns) {
        var l = data.length;
        if (rows === undefined) rows = l;
        if (columns === undefined) columns = rows;
        var min = Math.min(l, rows, columns);
        var matrix = Matrix.zeros(rows, columns);
        for (var i = 0; i < min; i++) {
            matrix[i][i] = data[i];
        }
        return matrix;
    }

    /**
     * Returns a matrix whose elements are the minimum between matrix1 and matrix2
     * @param matrix1
     * @param matrix2
     * @returns {Matrix}
     */
    static min(matrix1, matrix2) {
        var rows = matrix1.length;
        var columns = matrix1[0].length;
        var result = new Matrix(rows, columns);
        for (var i = 0; i < rows; i++) {
            for(var j = 0; j < columns; j++) {
                result[i][j] = Math.min(matrix1[i][j], matrix2[i][j]);
            }
        }
        return result;
    }

    /**
     * Returns a matrix whose elements are the maximum between matrix1 and matrix2
     * @param matrix1
     * @param matrix2
     * @returns {Matrix}
     */
    static max(matrix1, matrix2) {
        var rows = matrix1.length;
        var columns = matrix1[0].length;
        var result = new Matrix(rows, columns);
        for (var i = 0; i < rows; i++) {
            for(var j = 0; j < columns; j++) {
                result[i][j] = Math.max(matrix1[i][j], matrix2[i][j]);
            }
        }
        return result;
    }

    /**
     * Check that the provided value is a Matrix and tries to instantiate one if not
     * @param value - The value to check
     * @returns {Matrix}
     */
    static checkMatrix(value) {
        return Matrix.isMatrix(value) ? value : new Matrix(value);
    }

    /**
     * Returns true if the argument is a Matrix, false otherwise
     * @param value - The value to check
     * @return {boolean}
     */
    static isMatrix(value) {
        return (value != null) && (value.klass === 'Matrix');
    }

    /**
     * @property {number} - The number of elements in the matrix.
     */
    get size() {
        return this.rows * this.columns;
    }

    /**
     * Applies a callback for each element of the matrix. The function is called in the matrix (this) context.
     * @param {function} callback - Function that will be called with two parameters : i (row) and j (column)
     * @returns {Matrix} this
     */
    apply(callback) {
        if (typeof callback !== 'function') {
            throw new TypeError('callback must be a function');
        }
        var ii = this.rows;
        var jj = this.columns;
        for (var i = 0; i < ii; i++) {
            for (var j = 0; j < jj; j++) {
                callback.call(this, i, j);
            }
        }
        return this;
    }

    /**
     * Creates an exact and independent copy of the matrix
     * @returns {Matrix}
     */
    clone() {
        var newMatrix = new Matrix(this.rows, this.columns);
        for (var row = 0; row < this.rows; row++) {
            for (var column = 0; column < this.columns; column++) {
                newMatrix[row][column] = this[row][column];
            }
        }
        return newMatrix;
    }

    /**
     * Returns a new 1D array filled row by row with the matrix values
     * @returns {Array}
     */
    to1DArray() {
        var array = new Array(this.size);
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.columns; j++) {
                array[i * this.columns + j] = this[i][j];
            }
        }
        return array;
    }

    /**
     * Returns a 2D array containing a copy of the data
     * @returns {Array}
     */
    to2DArray() {
        var copy = new Array(this.rows);
        for (var i = 0; i < this.rows; i++) {
            copy[i] = [].concat(this[i]);
        }
        return copy;
    }

    /**
     * @returns {boolean} true if the matrix has one row
     */
    isRowVector() {
        return this.rows === 1;
    }

    /**
     * @returns {boolean} true if the matrix has one column
     */
    isColumnVector() {
        return this.columns === 1;
    }

    /**
     * @returns {boolean} true if the matrix has one row or one column
     */
    isVector() {
        return (this.rows === 1) || (this.columns === 1);
    }

    /**
     * @returns {boolean} true if the matrix has the same number of rows and columns
     */
    isSquare() {
        return this.rows === this.columns;
    }

    /**
     * @returns {boolean} true if the matrix is square and has the same values on both sides of the diagonal
     */
    isSymmetric() {
        if (this.isSquare()) {
            for (var i = 0; i < this.rows; i++) {
                for (var j = 0; j <= i; j++) {
                    if (this[i][j] !== this[j][i]) {
                        return false;
                    }
                }
            }
            return true;
        }
        return false;
    }

    /**
     * Sets a given element of the matrix. mat.set(3,4,1) is equivalent to mat[3][4]=1
     * @param {number} rowIndex - Index of the row
     * @param {number} columnIndex - Index of the column
     * @param {number} value - The new value for the element
     * @returns {Matrix} this
     */
    set(rowIndex, columnIndex, value) {
        this[rowIndex][columnIndex] = value;
        return this;
    }

    /**
     * Returns the given element of the matrix. mat.get(3,4) is equivalent to matrix[3][4]
     * @param {number} rowIndex - Index of the row
     * @param {number} columnIndex - Index of the column
     * @returns {number}
     */
    get(rowIndex, columnIndex) {
        return this[rowIndex][columnIndex];
    }

    /**
     * Fills the matrix with a given value. All elements will be set to this value.
     * @param {number} value - New value
     * @returns {Matrix} this
     */
    fill(value) {
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.columns; j++) {
                this[i][j] = value;
            }
        }
        return this;
    }

    /**
     * Negates the matrix. All elements will be multiplied by (-1)
     * @returns {Matrix} this
     */
    neg() {
        return this.mulS(-1);
    }

    /**
     * Returns a new array from the given row index
     * @param {number} index - Row index
     * @returns {Array}
     */
    getRow(index) {
        checkRowIndex(this, index);
        return [].concat(this[index]);
    }

    /**
     * Returns a new row vector from the given row index
     * @param {number} index - Row index
     * @returns {Matrix}
     */
    getRowVector(index) {
        return Matrix.rowVector(this.getRow(index));
    }

    /**
     * Sets a row at the given index
     * @param {number} index - Row index
     * @param {Array|Matrix} array - Array or vector
     * @returns {Matrix} this
     */
    setRow(index, array) {
        checkRowIndex(this, index);
        array = checkRowVector(this, array, true);
        this[index] = array;
        return this;
    }

    /**
     * Removes a row from the given index
     * @param {number} index - Row index
     * @returns {Matrix} this
     */
    removeRow(index) {
        checkRowIndex(this, index);
        if (this.rows === 1)
            throw new RangeError('A matrix cannot have less than one row');
        this.splice(index, 1);
        this.rows -= 1;
        return this;
    }

    /**
     * Adds a row at the given index
     * @param {number} [index = this.rows] - Row index
     * @param {Array|Matrix} array - Array or vector
     * @returns {Matrix} this
     */
    addRow(index, array) {
        if (array === undefined) {
            array = index;
            index = this.rows;
        }
        checkRowIndex(this, index, true);
        array = checkRowVector(this, array, true);
        this.splice(index, 0, array);
        this.rows += 1;
        return this;
    }

    /**
     * Swaps two rows
     * @param {number} row1 - First row index
     * @param {number} row2 - Second row index
     * @returns {Matrix} this
     */
    swapRows(row1, row2) {
        checkRowIndex(this, row1);
        checkRowIndex(this, row2);
        var temp = this[row1];
        this[row1] = this[row2];
        this[row2] = temp;
        return this;
    }

    /**
     * Returns a new array from the given column index
     * @param {number} index - Column index
     * @returns {Array}
     */
    getColumn(index) {
        checkColumnIndex(this, index);
        var column = new Array(this.rows);
        for (var i = 0; i < this.rows; i++) {
            column[i] = this[i][index];
        }
        return column;
    }

    /**
     * Returns a new column vector from the given column index
     * @param {number} index - Column index
     * @returns {Matrix}
     */
    getColumnVector(index) {
        return Matrix.columnVector(this.getColumn(index));
    }

    /**
     * Sets a column at the given index
     * @param {number} index - Column index
     * @param {Array|Matrix} array - Array or vector
     * @returns {Matrix} this
     */
    setColumn(index, array) {
        checkColumnIndex(this, index);
        array = checkColumnVector(this, array);
        for (var i = 0; i < this.rows; i++) {
            this[i][index] = array[i];
        }
        return this;
    }

    /**
     * Removes a column from the given index
     * @param {number} index - Column index
     * @returns {Matrix} this
     */
    removeColumn(index) {
        checkColumnIndex(this, index);
        if (this.columns === 1)
            throw new RangeError('A matrix cannot have less than one column');
        for (var i = 0; i < this.rows; i++) {
            this[i].splice(index, 1);
        }
        this.columns -= 1;
        return this;
    }

    /**
     * Adds a column at the given index
     * @param {number} [index = this.columns] - Column index
     * @param {Array|Matrix} array - Array or vector
     * @returns {Matrix} this
     */
    addColumn(index, array) {
        if (typeof array === 'undefined') {
            array = index;
            index = this.columns;
        }
        checkColumnIndex(this, index, true);
        array = checkColumnVector(this, array);
        for (var i = 0; i < this.rows; i++) {
            this[i].splice(index, 0, array[i]);
        }
        this.columns += 1;
        return this;
    }

    /**
     * Swaps two columns
     * @param {number} column1 - First column index
     * @param {number} column2 - Second column index
     * @returns {Matrix} this
     */
    swapColumns(column1, column2) {
        checkColumnIndex(this, column1);
        checkColumnIndex(this, column2);
        var temp, row;
        for (var i = 0; i < this.rows; i++) {
            row = this[i];
            temp = row[column1];
            row[column1] = row[column2];
            row[column2] = temp;
        }
        return this;
    }

    /**
     * Adds the values of a vector to each row
     * @param {Array|Matrix} vector - Array or vector
     * @returns {Matrix} this
     */
    addRowVector(vector) {
        vector = checkRowVector(this, vector);
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.columns; j++) {
                this[i][j] += vector[j];
            }
        }
        return this;
    }

    /**
     * Subtracts the values of a vector from each row
     * @param {Array|Matrix} vector - Array or vector
     * @returns {Matrix} this
     */
    subRowVector(vector) {
        vector = checkRowVector(this, vector);
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.columns; j++) {
                this[i][j] -= vector[j];
            }
        }
        return this;
    }

    /**
     * Multiplies the values of a vector with each row
     * @param {Array|Matrix} vector - Array or vector
     * @returns {Matrix} this
     */
    mulRowVector(vector) {
        vector = checkRowVector(this, vector);
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.columns; j++) {
                this[i][j] *= vector[j];
            }
        }
        return this;
    }

    /**
     * Divides the values of each row by those of a vector
     * @param {Array|Matrix} vector - Array or vector
     * @returns {Matrix} this
     */
    divRowVector(vector) {
        vector = checkRowVector(this, vector);
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.columns; j++) {
                this[i][j] /= vector[j];
            }
        }
        return this;
    }

    /**
     * Adds the values of a vector to each column
     * @param {Array|Matrix} vector - Array or vector
     * @returns {Matrix} this
     */
    addColumnVector(vector) {
        vector = checkColumnVector(this, vector);
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.columns; j++) {
                this[i][j] += vector[i];
            }
        }
        return this;
    }

    /**
     * Subtracts the values of a vector from each column
     * @param {Array|Matrix} vector - Array or vector
     * @returns {Matrix} this
     */
    subColumnVector(vector) {
        vector = checkColumnVector(this, vector);
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.columns; j++) {
                this[i][j] -= vector[i];
            }
        }
        return this;
    }

    /**
     * Multiplies the values of a vector with each column
     * @param {Array|Matrix} vector - Array or vector
     * @returns {Matrix} this
     */
    mulColumnVector(vector) {
        vector = checkColumnVector(this, vector);
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.columns; j++) {
                this[i][j] *= vector[i];
            }
        }
        return this;
    }

    /**
     * Divides the values of each column by those of a vector
     * @param {Array|Matrix} vector - Array or vector
     * @returns {Matrix} this
     */
    divColumnVector(vector) {
        vector = checkColumnVector(this, vector);
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.columns; j++) {
                this[i][j] /= vector[i];
            }
        }
        return this;
    }

    /**
     * Multiplies the values of a row with a scalar
     * @param {number} index - Row index
     * @param {number} value
     * @returns {Matrix} this
     */
    mulRow(index, value) {
        checkRowIndex(this, index);
        for (var i = 0; i < this.columns; i++) {
            this[index][i] *= value;
        }
        return this;
    }

    /**
     * Multiplies the values of a column with a scalar
     * @param {number} index - Column index
     * @param {number} value
     * @returns {Matrix} this
     */
    mulColumn(index, value) {
        checkColumnIndex(this, index);
        for (var i = 0; i < this.rows; i++) {
            this[i][index] *= value;
        }
    }

    /**
     * Returns the maximum value of the matrix
     * @returns {number}
     */
    max() {
        var v = this[0][0];
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.columns; j++) {
                if (this[i][j] > v) {
                    v = this[i][j];
                }
            }
        }
        return v;
    }

    /**
     * Returns the index of the maximum value
     * @returns {Array}
     */
    maxIndex() {
        var v = this[0][0];
        var idx = [0, 0];
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.columns; j++) {
                if (this[i][j] > v) {
                    v = this[i][j];
                    idx[0] = i;
                    idx[1] = j;
                }
            }
        }
        return idx;
    }

    /**
     * Returns the minimum value of the matrix
     * @returns {number}
     */
    min() {
        var v = this[0][0];
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.columns; j++) {
                if (this[i][j] < v) {
                    v = this[i][j];
                }
            }
        }
        return v;
    }

    /**
     * Returns the index of the minimum value
     * @returns {Array}
     */
    minIndex() {
        var v = this[0][0];
        var idx = [0, 0];
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.columns; j++) {
                if (this[i][j] < v) {
                    v = this[i][j];
                    idx[0] = i;
                    idx[1] = j;
                }
            }
        }
        return idx;
    }

    /**
     * Returns the maximum value of one row
     * @param {number} row - Row index
     * @returns {number}
     */
    maxRow(row) {
        checkRowIndex(this, row);
        var v = this[row][0];
        for (var i = 1; i < this.columns; i++) {
            if (this[row][i] > v) {
                v = this[row][i];
            }
        }
        return v;
    }

    /**
     * Returns the index of the maximum value of one row
     * @param {number} row - Row index
     * @returns {Array}
     */
    maxRowIndex(row) {
        checkRowIndex(this, row);
        var v = this[row][0];
        var idx = [row, 0];
        for (var i = 1; i < this.columns; i++) {
            if (this[row][i] > v) {
                v = this[row][i];
                idx[1] = i;
            }
        }
        return idx;
    }

    /**
     * Returns the minimum value of one row
     * @param {number} row - Row index
     * @returns {number}
     */
    minRow(row) {
        checkRowIndex(this, row);
        var v = this[row][0];
        for (var i = 1; i < this.columns; i++) {
            if (this[row][i] < v) {
                v = this[row][i];
            }
        }
        return v;
    }

    /**
     * Returns the index of the maximum value of one row
     * @param {number} row - Row index
     * @returns {Array}
     */
    minRowIndex(row) {
        checkRowIndex(this, row);
        var v = this[row][0];
        var idx = [row, 0];
        for (var i = 1; i < this.columns; i++) {
            if (this[row][i] < v) {
                v = this[row][i];
                idx[1] = i;
            }
        }
        return idx;
    }

    /**
     * Returns the maximum value of one column
     * @param {number} column - Column index
     * @returns {number}
     */
    maxColumn(column) {
        checkColumnIndex(this, column);
        var v = this[0][column];
        for (var i = 1; i < this.rows; i++) {
            if (this[i][column] > v) {
                v = this[i][column];
            }
        }
        return v;
    }

    /**
     * Returns the index of the maximum value of one column
     * @param {number} column - Column index
     * @returns {Array}
     */
    maxColumnIndex(column) {
        checkColumnIndex(this, column);
        var v = this[0][column];
        var idx = [0, column];
        for (var i = 1; i < this.rows; i++) {
            if (this[i][column] > v) {
                v = this[i][column];
                idx[0] = i;
            }
        }
        return idx;
    }

    /**
     * Returns the minimum value of one column
     * @param {number} column - Column index
     * @returns {number}
     */
    minColumn(column) {
        checkColumnIndex(this, column);
        var v = this[0][column];
        for (var i = 1; i < this.rows; i++) {
            if (this[i][column] < v) {
                v = this[i][column];
            }
        }
        return v;
    }

    /**
     * Returns the index of the minimum value of one column
     * @param {number} column - Column index
     * @returns {Array}
     */
    minColumnIndex(column) {
        checkColumnIndex(this, column);
        var v = this[0][column];
        var idx = [0, column];
        for (var i = 1; i < this.rows; i++) {
            if (this[i][column] < v) {
                v = this[i][column];
                idx[0] = i;
            }
        }
        return idx;
    }

    /**
     * Returns an array containing the diagonal values of the matrix
     * @returns {Array}
     */
    diag() {
        var min = Math.min(this.rows, this.columns);
        var diag = new Array(min);
        for (var i = 0; i < min; i++) {
            diag[i] = this[i][i];
        }
        return diag;
    }

    /**
     * Returns the sum of all elements of the matrix
     * @returns {number}
     */
    sum() {
        var v = 0;
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.columns; j++) {
                v += this[i][j];
            }
        }
        return v;
    }

    /**
     * Returns the mean of all elements of the matrix
     * @returns {number}
     */
    mean() {
        return this.sum() / this.size;
    }

    /**
     * Returns the product of all elements of the matrix
     * @returns {number}
     */
    prod() {
        var prod = 1;
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.columns; j++) {
                prod *= this[i][j];
            }
        }
        return prod;
    }

    /**
     * Computes the cumulative sum of the matrix elements (in place, row by row)
     * @returns {Matrix} this
     */
    cumulativeSum() {
        var sum = 0;
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.columns; j++) {
                sum += this[i][j];
                this[i][j] = sum;
            }
        }
        return this;
    }

    /**
     * Computes the dot (scalar) product between the matrix and another
     * @param {Matrix} vector2 vector
     * @returns {number}
     */
    dot(vector2) {
        if (Matrix.isMatrix(vector2)) vector2 = vector2.to1DArray();
        var vector1 = this.to1DArray();
        if (vector1.length !== vector2.length) {
            throw new RangeError('vectors do not have the same size');
        }
        var dot = 0;
        for (var i = 0; i < vector1.length; i++) {
            dot += vector1[i] * vector2[i];
        }
        return dot;
    }

    /**
     * Returns the matrix product between this and other
     * @returns {Matrix}
     */
    mmul(other) {
        other = Matrix.checkMatrix(other);
        if (this.columns !== other.rows)
            console.warn('Number of columns of left matrix are not equal to number of rows of right matrix.');

        var m = this.rows;
        var n = this.columns;
        var p = other.columns;

        var result = new Matrix(m, p);

        var Bcolj = new Array(n);
        for (var j = 0; j < p; j++) {
            for (var k = 0; k < n; k++)
                Bcolj[k] = other[k][j];

            for (var i = 0; i < m; i++) {
                var Arowi = this[i];

                var s = 0;
                for (k = 0; k < n; k++)
                    s += Arowi[k] * Bcolj[k];

                result[i][j] = s;
            }
        }
        return result;
    }

    /**
     * Transposes the matrix and returns a new one containing the result
     * @returns {Matrix}
     */
    transpose() {
        var result = new Matrix(this.columns, this.rows);
        for (var i = 0; i < this.rows; i++) {
            for (var j = 0; j < this.columns; j++) {
                result[j][i] = this[i][j];
            }
        }
        return result;
    }

    /**
     * Sorts the rows (in place)
     * @param {function} compareFunction - usual Array.prototype.sort comparison function
     * @returns {Matrix} this
     */
    sortRows(compareFunction) {
        if (compareFunction === undefined) compareFunction = compareNumbers;
        for (var i = 0; i < this.rows; i++) {
            this[i].sort(compareFunction);
        }
        return this;
    }

    /**
     * Sorts the columns (in place)
     * @param {function} compareFunction - usual Array.prototype.sort comparison function
     * @returns {Matrix} this
     */
    sortColumns(compareFunction) {
        if (compareFunction === undefined) compareFunction = compareNumbers;
        for (var i = 0; i < this.columns; i++) {
            this.setColumn(i, this.getColumn(i).sort(compareFunction));
        }
        return this;
    }

    /**
     * Returns a subset of the matrix
     * @param {number} startRow - First row index
     * @param {number} endRow - Last row index
     * @param {number} startColumn - First column index
     * @param {number} endColumn - Last column index
     * @returns {Matrix}
     */
    subMatrix(startRow, endRow, startColumn, endColumn) {
        if ((startRow > endRow) || (startColumn > endColumn) || (startRow < 0) || (startRow >= this.rows) || (endRow < 0) || (endRow >= this.rows) || (startColumn < 0) || (startColumn >= this.columns) || (endColumn < 0) || (endColumn >= this.columns)) {
            throw new RangeError('Argument out of range');
        }
        var newMatrix = new Matrix(endRow - startRow + 1, endColumn - startColumn + 1);
        for (var i = startRow; i <= endRow; i++) {
            for (var j = startColumn; j <= endColumn; j++) {
                newMatrix[i - startRow][j - startColumn] = this[i][j];
            }
        }
        return newMatrix;
    }

    /**
     * Returns a subset of the matrix based on an array of row indices
     * @param {Array} indices - Array containing the row indices
     * @param {number} [startColumn = 0] - First column index
     * @param {number} [endColumn = this.columns-1] - Last column index
     * @returns {Matrix}
     */
    subMatrixRow(indices, startColumn, endColumn) {
        if (startColumn === undefined) startColumn = 0;
        if (endColumn === undefined) endColumn = this.columns - 1;
        if ((startColumn > endColumn) || (startColumn < 0) || (startColumn >= this.columns) || (endColumn < 0) || (endColumn >= this.columns)) {
            throw new RangeError('Argument out of range');
        }

        var newMatrix = new Matrix(indices.length, endColumn - startColumn + 1);
        for (var i = 0; i < indices.length; i++) {
            for (var j = startColumn; j <= endColumn; j++) {
                if (indices[i] < 0 || indices[i] >= this.rows) {
                    throw new RangeError('Row index out of range: ' + indices[i]);
                }
                newMatrix[i][j - startColumn] = this[indices[i]][j];
            }
        }
        return newMatrix;
    }

    /**
     * Returns a subset of the matrix based on an array of column indices
     * @param {Array} indices - Array containing the column indices
     * @param {number} [startRow = 0] - First row index
     * @param {number} [endRow = this.rows-1] - Last row index
     * @returns {Matrix}
     */
    subMatrixColumn(indices, startRow, endRow) {
        if (startRow === undefined) startRow = 0;
        if (endRow === undefined) endRow = this.rows - 1;
        if ((startRow > endRow) || (startRow < 0) || (startRow >= this.rows) || (endRow < 0) || (endRow >= this.rows)) {
            throw new RangeError('Argument out of range');
        }

        var newMatrix = new Matrix(endRow - startRow + 1, indices.length);
        for (var i = 0; i < indices.length; i++) {
            for (var j = startRow; j <= endRow; j++) {
                if (indices[i] < 0 || indices[i] >= this.columns) {
                    throw new RangeError('Column index out of range: ' + indices[i]);
                }
                newMatrix[j - startRow][i] = this[j][indices[i]];
            }
        }
        return newMatrix;
    }

    /**
     * Returns the trace of the matrix (sum of the diagonal elements)
     * @returns {number}
     */
    trace() {
        var min = Math.min(this.rows, this.columns);
        var trace = 0;
        for (var i = 0; i < min; i++) {
            trace += this[i][i];
        }
        return trace;
    }
}

Matrix.prototype.klass = 'Matrix';

module.exports = Matrix;

/**
 * @private
 * Check that a row index is not out of bounds
 * @param {Matrix} matrix
 * @param {number} index
 * @param {boolean} [outer]
 */
function checkRowIndex(matrix, index, outer) {
    var max = outer ? matrix.rows : matrix.rows - 1;
    if (index < 0 || index > max)
        throw new RangeError('Row index out of range');
}

/**
 * @private
 * Check that the provided vector is an array with the right length
 * @param {Matrix} matrix
 * @param {Array|Matrix} vector
 * @param {boolean} copy
 * @returns {Array}
 * @throws {RangeError}
 */
function checkRowVector(matrix, vector, copy) {
    if (Matrix.isMatrix(vector)) {
        vector = vector.to1DArray();
    } else if (copy) {
        vector = [].concat(vector);
    }
    if (vector.length !== matrix.columns)
        throw new RangeError('vector size must be the same as the number of columns');
    return vector;
}

/**
 * @private
 * Check that the provided vector is an array with the right length
 * @param {Matrix} matrix
 * @param {Array|Matrix} vector
 * @param {boolean} copy
 * @returns {Array}
 * @throws {RangeError}
 */
function checkColumnVector(matrix, vector, copy) {
    if (Matrix.isMatrix(vector)) {
        vector = vector.to1DArray();
    } else if (copy) {
        vector = [].concat(vector);
    }
    if (vector.length !== matrix.rows)
        throw new RangeError('vector size must be the same as the number of rows');
    return vector;
}

/**
 * @private
 * Check that a column index is not out of bounds
 * @param {Matrix} matrix
 * @param {number} index
 * @param {boolean} [outer]
 */
function checkColumnIndex(matrix, index, outer) {
    var max = outer ? matrix.columns : matrix.columns - 1;
    if (index < 0 || index > max)
        throw new RangeError('Column index out of range');
}

/**
 * @private
 * Check that two matrices have the same dimensions
 * @param {Matrix} matrix
 * @param {Matrix} otherMatrix
 */
function checkDimensions(matrix, otherMatrix) {
    if (matrix.rows !== otherMatrix.length ||
        matrix.columns !== otherMatrix[0].length) {
        throw new RangeError('Matrices dimensions must be equal');
    }
}

function compareNumbers(a, b) {
    return a - b;
}

/*
Synonyms
 */

Matrix.random = Matrix.rand;
Matrix.diagonal = Matrix.diag;
Matrix.prototype.diagonal = Matrix.prototype.diag;
Matrix.identity = Matrix.eye;
Matrix.prototype.negate = Matrix.prototype.neg;

/*
Add dynamically instance and static methods for mathematical operations
 */

var inplaceOperator = `
(function %name%(value) {
    if (typeof value === 'number') return this.%name%S(value);
    return this.%name%M(value);
})
`;

var inplaceOperatorScalar = `
(function %name%S(value) {
    for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
            this[i][j] = this[i][j] %op% value;
        }
    }
    return this;
})
`;

var inplaceOperatorMatrix = `
(function %name%M(matrix) {
    checkDimensions(this, matrix);
    for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
            this[i][j] = this[i][j] %op% matrix[i][j];
        }
    }
    return this;
})
`;

var staticOperator = `
(function %name%(matrix, value) {
    var newMatrix = new Matrix(matrix);
    return newMatrix.%name%(value);
})
`;

var inplaceMethod = `
(function %name%() {
    for (var i = 0; i < this.rows; i++) {
        for (var j = 0; j < this.columns; j++) {
            this[i][j] = %method%(this[i][j]);
        }
    }
    return this;
})
`;

var staticMethod = `
(function %name%(matrix) {
    var newMatrix = new Matrix(matrix);
    return newMatrix.%name%();
})
`;

var operators = [
    // Arithmetic operators
    ['+', 'add'],
    ['-', 'sub', 'subtract'],
    ['*', 'mul', 'multiply'],
    ['/', 'div', 'divide'],
    ['%', 'mod', 'modulus'],
    // Bitwise operators
    ['&', 'and'],
    ['|', 'or'],
    ['^', 'xor'],
    ['<<', 'leftShift'],
    ['>>', 'signPropagatingRightShift'],
    ['>>>', 'rightShift', 'zeroFillRightShift']
];

for (var operator of operators) {
    for (var i = 1; i < operator.length; i++) {
        Matrix.prototype[operator[i]] = eval(fillTemplateFunction(inplaceOperator, {name: operator[i], op: operator[0]}));
        Matrix.prototype[operator[i] + 'S'] = eval(fillTemplateFunction(inplaceOperatorScalar, {name: operator[i] + 'S', op: operator[0]}));
        Matrix.prototype[operator[i] + 'M'] = eval(fillTemplateFunction(inplaceOperatorMatrix, {name: operator[i] + 'M', op: operator[0]}));

        Matrix[operator[i]] = eval(fillTemplateFunction(staticOperator, {name: operator[i]}));
    }
}

var methods = [
    ['~', 'not']
];

[
    'abs', 'acos', 'acosh', 'asin', 'asinh', 'atan', 'atanh', 'cbrt', 'ceil',
    'clz32', 'cos', 'cosh', 'exp', 'expm1', 'floor', 'fround', 'log', 'log1p',
    'log10', 'log2', 'round', 'sign', 'sin', 'sinh', 'sqrt', 'tan', 'tanh', 'trunc'
].forEach(function (mathMethod) {
    methods.push(['Math.' + mathMethod, mathMethod]);
});

for (var method of methods) {
    for (var i = 1; i < method.length; i++) {
        Matrix.prototype[method[i]] = eval(fillTemplateFunction(inplaceMethod, {name: method[i], method: method[0]}));
        Matrix[method[i]] = eval(fillTemplateFunction(staticMethod, {name: method[i]}));
    }
}

function fillTemplateFunction(template, values) {
    for (var i in values) {
        template = template.replace(new RegExp('%' + i + '%', 'g'), values[i]);
    }
    return template;
}

},{}]},{},[3]);
