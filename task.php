<!doctype html>
<html>

<?php
  session_start();

  include($_SERVER["DOCUMENT_ROOT"]."/code/php/AC.php");
  $user_name = check_logged(); /// function checks if visitor is logged in.
  $admin = false;

  if ($user_name == "") {
    // user is not logged in

  } else {
    $admin = true;
    echo('<script type="text/javascript"> user_name = "'.$user_name.'"; </script>'."\n");
    echo('<script type="text/javascript"> admin = '.($admin?"true":"false").'; </script>'."\n");
  }

  $subjid = "";
  $sessionid = "";
  if( isset($_SESSION['ABCD']) && isset($_SESSION['ABCD']['delay-discounting']) ) {
     if (isset($_SESSION['ABCD']['delay-discounting']['subjid'])) {
        $subjid  = $_SESSION['ABCD']['delay-discounting']['subjid'];
     }
     if (isset($_SESSION['ABCD']['delay-discounting']['sessionid'])) {
        $sessionid  = $_SESSION['ABCD']['delay-discounting']['sessionid'];
     }
  }
  echo('<script type="text/javascript"> SubjectID = "'.$subjid.'"; </script>'."\n");
  echo('<script type="text/javascript"> Session = "'.$sessionid.'"; </script>'."\n");

   $permissions = list_permissions_for_user( $user_name );

   $site = "";
   foreach ($permissions as $per) {
     $a = explode("Site", $per); // permissions should be structured as "Site<site name>"

     if (count($a) > 0) {
        $site = $a[1];
	break;
     }
   }
   if ($site == "") {
     echo (json_encode ( array( "message" => "Error: no site assigned to this user" ) ) );
     return;
   }
   echo('<script type="text/javascript"> Site = "'.$site.'"; </script>'."\n");

?>


  <head>
    <title>Delay-Discounting Task</title>
    <meta charset="utf-8" />
    <!-- Load jQuery -->
    <script src="js/jquery.min.js"></script>
    <!-- <script src="/js/jquery.ui.touch-punch.min.js"></script>     -->
    <script src='js/moment.min.js'></script>
   
    <!-- Load the jspsych library and plugins -->
    <script src="js/jspsych/jspsych.js"></script>
    <script src="js/jspsych/plugins/jspsych-text.js"></script>
    <script src="js/jspsych/plugins/jspsych-single-stim.js"></script>
    <script src="js/jspsych/plugins/jspsych-button-response.js"></script>
    <!-- Load plotting library -->
    <script src="js/plotly-latest.min.js"></script>
    
    <!-- Load the stylesheet -->
    <!-- <link href="experiment.css" type="text/css" rel="stylesheet"></link> -->
    <link href="js/jspsych/css/jspsych.css" rel="stylesheet" type="text/css"></link>
    <link href='https://fonts.googleapis.com/css?family=Lato:300' rel='stylesheet' type='text/css'>
    <link href='https://fonts.googleapis.com/css?family=Open+Sans:300italic' rel='stylesheet' type='text/css'>
<style>
body {
  backgroud-color: black;
  color: white;
}
h1 {
   color: #ffffff;
   font-family: 'Lato', sans-serif;
   font-size: 54px;
   font-weight: 300;
   line-height: 58px;
   margin: 0 0 58px;
   border-bottom: double #555;
   padding-bottom: 30px;
}
h2 {
   color: #ffffff;
   font-family: 'Lato', sans-serif;
   font-size: 34px;
   font-weight: 300;
   line-height: 48px;
   margin: 0 0 48px;
   padding-bottom: 30px;
}
p {
   color: #adb7bd;
   font-family: 'Open Sans', Arial, sans-serif;
   font-size: 16px;
   line-height: 26px;
   text-indent: 0px;
    margin: 0;
    font-size: 32px;
    line-height: 1.2em;
}
a {
   color: #fe921f;
   text-decoration: underline;
}
a:hover { color: #ffffff }
.date {
      background: #fe921f;
      color: #ffffff;
      display: inline-block;
      font-family: 'Lato', sans-serif;
      font-size: 12px;
      font-weight: bold;
      line-height: 12px;
      letter-spacing: 1px;
      margin: 0 0 30px;
      padding: 10px 15px 8px;
      text-transform: uppercase;
}

.date2 { color: #bbc3c8; background: #292929; display: inline-block; font-family: 'Georgia', serif; font-style: italic; font-size: 18px; line-height: 22px; margin: 0 0 20px 18px; padding: 10px 12px 8px; position: absolute; bottom: -36px; }

.jspsych-btn {
  position: absolute;
  bottom: 40px;
  border-radius: 40px;
  width: 80px;
  height: 80px;
  font-size: 32pt;
  color: gray;
  
}
#jspsych-button-response-button-0 {
  box-shadow: 0px 0px 8px #fff;
  left: 30%;
}
#jspsych-button-response-button-1 {
  box-shadow: 0px 0px 8px #fff;
  right: 30%;


}

    .left {
	position: absolute;
	top: 200px;
	left: 20%;
	font-size: 34px;
    }
    .right {
	position: absolute;
	top: 200px;
	left: 60%;
	font-size: 34px;
    }
    .crosshair {
	font-size: 60px;
	position: absolute;
	left: 50%;
	top: 200px;
    }
</style>


  </head>

  <body bgcolor="#292929">
    <div id="jspsych_target"></div>
  </body>
  
  <!-- Load the math libraries to calculate scores -->
  <script src="js/ml-matrix-bundle.js"></script>
  <script src="js/curve-fit-bundle.js"></script>

  <script>

function exportToCsv(filename, rows) {
    var k = { "SubjectID": 1, "Site": 1, "Session": 1 };
    for (var i = 0; i < rows.length; i++) {
       var k2 = Object.keys(rows[i]);
       for (var j = 0; j < k2.length; j++) {
          k[k2[j]] = 1;
       } 
    }
    k = Object.keys(k);

    var csvFile = k.join(",") + "\n";
    for (var i = 0; i < rows.length; i++) {
       rows[i]['SubjectID'] = SubjectID;
       rows[i]['Site'] = Site;
       rows[i]['Session'] = Session;
       csvFile += k.map(function(a) { return rows[i][a] }).join(",") + "\n";
    }
    
    var blob = new Blob([csvFile], { type: 'text/csv;charset=utf-8;' });
    if (navigator.msSaveBlob) { // IE 10+
	navigator.msSaveBlob(blob, filename);
    } else {
	var link = document.createElement("a");
	if (link.download !== undefined) { // feature detection
	    // Browsers that support HTML5 download attribute
	    var url = URL.createObjectURL(blob);
	    link.setAttribute("href", url);
	    link.setAttribute("download", filename);
	    link.style.visibility = 'hidden';
	    document.body.appendChild(link);
	    link.click();
	    document.body.removeChild(link);
	}
    }
}



    var post_trial_gap = function() {
        return Math.floor( Math.random() * 1000 ) + 500;
    }


    // we need to query for key "1"/"left cursor" or "6"/"right cursor" 
    var instructions = "<p>In this game, you will be asked to make some choices between getting some amount of money right now or waiting to get a larger amount of money in the future.</p><br/><p>These money amounts are pretend. You won't actually get these amounts of money for this game, but we ask you to choose between the amounts of money as if they were real.</p><br/><p>You need to use the L and R buttons to make each choice, and you can take as much time as you need to make your choice.</p><br/><p>Press the button to begin.</p>";

    var instructions2 = "<p>Press '8' for left choice, press '9' for right choice.</p>";

    var thanks = "<center><p>Thank you for participating!</p></center>";

//
// The following code has been translated from python code originally written by
// Micky Koffarnus at Virginia Tech.
//
// settings (Discounting-ABCD($100).py, AdjAmt discounting everything.py)
sid = SubjectID;
exp = Session;
probablity = 0;
losses = 0;
past = 0;
explicit0 = 0;
commodityD = "$";
commodityI = "$";
amountD = "100";
amountI = "100";
directory = "data/" + Site;
customdelays = 'n';
// Original delays contained a 25/9131.25 years option, possible would also be a 6 hour question			
// delays=['1 day', '1 week', '1 month', '3 months', '1 year', '5 years', '25 years'];
// xx = [1,7,30.44,91.32,365.25,1826.25,9131.25];
delays=['1 day', '1 week', '1 month', '3 months', '1 year', '5 years'];
xx = [1,7,30.44,91.32,365.25,1826.25];
tr_len = 2; // how long each stimulus is displayed?
x1 = Array.apply(null, { length: 5 }).map(Number.call, Number).map(function(a) { return a+1; });
x2 = x1.map(function(a) { return a/2.0; });

// create 1168 random samples from x2
start = 0;
end = 1168;
isi_array = Array.apply(null, { length: end }).map(function(a) { return Math.floor(Math.random() * (x2.length - start + 1) + start); });
trs=1168*tr_len+ isi_array.reduce(function(a, b) { return a+b; }, 0) *2; //  up to 1168 questions
gn_sec_n=trs*tr_len; //  total time for presentation
gn_keystroke = 0;
question=0;
currenttask=0;
sorteddelays=[];
for (var i = 0; i < delays.length; i++) {
   sorteddelays.push( delays[i] );
}
currentdelay=0;
delays = jsPsych.randomization.shuffle( delays );  // shuffle(delays)
numdelays=delays.length;
qIndex=1;
screenText=['',''];
screenText2=['',''];
ip=0;
ii=0;
iv=0;
subQ=-1;
krow=0;
blocknum=0;
// log file name
//if not os.path.exists(('%s\\data') % (directory)):
//os.makedirs(('%s\\data') % (directory))
//if not os.path.exists('%s\\data\\%s' % (directory, sid)):
//os.makedirs('%s\\data\\%s' % (directory, sid))
//log_filename = '%s\\data\\%s\\AdjAmt_%s_%s_' % (directory, sid, sid, exp) + time.strftime ('%m-%d-%Y_%Hh-%Mm.csv')
//pickle_filename = '%s\\data\\%s\\AdjAmt_%s_%s.p' % (directory, sid, sid, exp)
isi_array = jsPsych.randomization.shuffle(isi_array); // shuffle(isi_array)
NEW_BLOCK_TIME=5;
isi_array.unshift(NEW_BLOCK_TIME);
isi = 0;
isi_index=0;
now=1;
firsttime=0;
response=0;
stim_onset=0;
fixon=0;
newseton=1;
taskendon=0;
trialnum=0;
goodrow = [];
IPs = Array.apply(null, { length: numdelays }).map(function(a) { return 0.0; }); // [0.0]*numdelays
amount=0.5;
fontsize = 60;
amtText=['',''];
amtText0I='';
amtText0D='';
var k=0;
var JBpass="FAIL";
var CRC01 = "FAIL";
var CRC02 = "FAIL";
var JB1=0;
var JB2=0;
		    
// get a string for this number that has a comma at the thousands place
function group(number) {
    return number.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
}

// this function will create the left/right choices for each iteration, uses all of these global variables
function getState(t) {

    //if (t > isi+isi_array[isi_index]) {
	newseton=0;
	fixon = 0;
	taskendon=0;
	if (firsttime) {
	    now = Math.round(Math.random()); // will be 0 or 1
	    stim_onset=t;
	    firsttime=0;
	}
	if (question < 6) {
	    delay = delays[currentdelay];
	    amtText0I = "$0";
	    if (amountI < 1000) {
		amtText[now] = "$" + (amountI*amount).toFixed(2);
	    } else {
		amtText[now] = "$" + group(Math.floor(amountI*amount));
	    }
	    amtText0D = "$0";
	    if (amountD<1000) {
		amtText[1-now] = "$" + (amountD);
	    } else {
		amtText[1-now] = "$" + group(Math.floor(amountD));
	    }
	    screenText[now] = amtText[now];
	    screenText[1-now] = amtText[1-now];
	    screenText2[now] = "now";
	    screenText2[1-now] = "in " + delay;
	    
            screenText[now] = "get " + screenText[now];
	    screenText[1-now] = "get " + screenText[1-now];
	    
            if (gn_keystroke > 0) {
		firsttime = 1;
		fixon = 1;
		if ((gn_keystroke == 1) && (now == 0)) {
		    response = 0;
		} else if((gn_keystroke == 3) && (now == 1)) {
		    response = 0;
		} else {
		    response = 1;
		}
		isi = t;
		isi_index = isi_index + 1;
		screenText[0] = "";
		screenText[1] = "";
		screenText2[0] = "";
		screenText2[1] = "";
		// logfile.write("%i,%f,%f,%f,%s,%i,%i\n" % (trialnum, stim_onset, t, amountI*amount, delays[currentdelay], response, now))
		console.log(trialnum + "," + stim_onset + "," + t + "," + amountI*amount + "," + delays[currentdelay] + "," + response + "," + now);
		
		trialnum = trialnum+1;
		if (response == losses) {
		    amount = amount - (Math.pow(0.5, question + 2));
		} else {
		    amount = amount + (Math.pow(0.5, question + 2));
		}
		if (question < 5) {
		    question = question + 1;
		} else {
		    console.log("amount: " + amount);
		    IPs[sorteddelays.indexOf(delays[currentdelay])] = amount;
		    console.log("IPs: " + IPs.join(","));
		    question = 0;
		    amount = 0.5;
		    gn_keystroke = 0;
		    if (currentdelay == numdelays-1) {
			console.log("sorted delays: " + sorteddelays.join(","));
			console.log("IPs: " + IPs.join(","));

			JB1 = 0;
			for (var i = 0; i < numdelays-1; i++) {
			    if (IPs[i+1]-IPs[i] > 0.2)
				JB1 = JB1 + 1;
			}
			JB2 = 0;
			if (IPs[0] - IPs[numdelays-1] < 0.1) {
			    JB2 = 1;
			}
			JBpass = "PASS";
			CRC01 = "PASS";
			CRC02 = "PASS";
			if (JB1 > 1) {
			       JBpass = "FAIL";
			       CRC01 = "FAIL";
			}
			if (JB2 > 0) {
			  JBpass = "FAIL";
			  CRC02 = "FAIL";
			}
			console.log("JB Rule 1, " + JB1);
			console.log("JB Rule 2, " + JB2);
			
                        //xvalues = numpy.array(xx)
			//yvalues = numpy.array(IPs)
			xvalues = xx;
			yvalues = IPs;
			var data = [];
			for (var i = 0; i < xx.length; i++) {
			   data.push([xvalues[i], yvalues[i]]);
		        }
			var f = function(x, k) {
			    var result = cf.getMatrix(x.rows, x.columns);
			    for (var i = 0; i < x.rows; i++)
				result[i][0] = 1/(1+(x[i][0]*k[0]));
			    return result;
			}
		     
			cf.curve_fit(data, f, [-2], [2]);
			console.log("Consistent responding check (0.1 criterion): " + JB1);
			console.log("Consistent responding check (0.2 criterion): " + JB2);
			console.log("Consistency: " + JBpass);
			console.log("k: " + cf.params[0][0]);
			k=cf.params[0][0];
			console.log("ln(k) : " + Math.log(cf.params[0][0]));
			//popt, pconv = curve_fit(func, xvalues, yvalues, p0 = 0.01)
			//screenText2[0] = "Consistency: %s" % (JBpass)
			//screenText[1] = "k value: %2.4f" % (float(popt))
			//screenText2[1] = "ln(k) value: %2.4f" % (float(log(popt)))
			//logfile.write("k value, %f\n" % float(popt))
			//logfile.write("ln(k) value, %f\n" % float(log(popt)))
			//IPs.append(popt)
			IPs.push(JB1);
			IPs.push(JB2);
			//pickle.dump(IPs, open(pickle_filename, "wb"))
			taskendon=1;
			// taskend.parameters.text = "Task Complete";
			console.log("Task Complete");
			fixon=0;
			isi=t+1000;
						    
		    } else {
			currentdelay=currentdelay+1;
			isi=t+NEW_BLOCK_TIME;
			newseton=1;
			// newset.parameters.text = "Next delay: " + (delays[currentdelay]);
			console.log("Next delay: " + delays[currentdelay]);
			fixon=0;
		    }
		}
	    }
	}
    //} else {
//	firsttime = 1;
  //  }
    gn_keystroke = 0;
    return 1;
}

// a single trial is a text stimulus
var trial = {
    type: 'button-response',
    is_html: true,
    timing_post_trial: 0,
    stimulus: function() {
	// here we need to call getState, but this is before we actually show the next one, we need to call this after we did this one
	gn_keystroke = 0;
	getState(t); // set the text strings
	return "<div class=\"left\">" + screenText[0] + " " + screenText2[0] + "</div> <div class=\"right\">" + screenText[1] + " " + screenText2[1] + "</div>";
    }
}
var crosshair = {
    type: 'single-stim',
    is_html: true,
    timing_post_trial: 0,
    choices: 'none',
    timing_stim: [400],
    timing_response: 400,
    stimulus: ["<div class=\"crosshair\">+</div>"]
}

// all trials are at most 1168
var all_trials = Array.apply(null, { length: 1168 }).map(function(a) { return [crosshair, trial]; }).reduce(function(a, b) { return a.concat(b); });

t = 6; // this is a time - but we control time differently here... we should not use that t parameter
var test_block = {
    type: 'single-stim',
    choices: ['L','R'],
    timeline: all_trials,
    on_finish: function(data) {
	if (data.key_press == jsPsych.pluginAPI.convertKeyCharacterToKeyCode('esc'))
	    jsPsych.endCurrentTimeline();
	jsPsych.data.addDataToLastTrial({
	    choice1: screenText[0] + " " + screenText2[0],
	    choice2: screenText[1] + " " + screenText2[1],
	    trialnum: trialnum,
	    stim_onset: stim_onset,
	    t: jsPsych.totalTime() / 1000.0,
	    amount: amountI*amount,
	    delays: delays[currentdelay],
	    response: response,
	    now: now
	});

	// The next stimulus needs to change, call getState and tell it what happened on the way here
	gn_keystroke = 0; // call once to set the next iteration values
	if (data.button_pressed == 0)
	    gn_keystroke = 1;
	if (data.button_pressed == 1)
	    gn_keystroke = 3;
	t = jsPsych.totalTime() / 1000;
	var a = getState(t); // compute new amounts
	if (taskendon == 1) {
	    jsPsych.endCurrentTimeline();
	}
    }
};

    var timeline = [];
    timeline.push( { type: 'button-response', choices: ['>'], is_html: true,
	     stimulus: instructions } );
    //timeline.push( { type: 'text', text: instructions2 } );
    timeline.push( test_block );
    timeline.push( { type: 'button-response', choices: ['->'], is_html: true, stimulus: thanks } ); 

    jsPsych.init({
	timeline: timeline,
	on_finish: function(data) {
	    // we should store the delays and IPs as well
	    var d = {};
	    sorteddelays.forEach(function( val, key) { d["IP "+val] = IPs[key]; });
	    d.k = k;
	    d.logk = Math.log(k);
	    d.Consistency = JBpass;
	    d.Consistent_resp_check01 = CRC01;
	    d.Consistent_resp_check02 = CRC02;
	    jsPsych.data.addDataToLastTrial(d);

	    ud = makeUnique( jsPsych.data.getData(), 'ded_' );
			  
	    jQuery.post('code/php/events.php',
			{ "data": JSON.stringify(ud), "date": moment().format() }, function(data) {
			    if (typeof data.ok == 'undefined' || data.ok == 0) {
				//  alert('Error: ' + data.message);
			    }
			    // export as csv for download on client
			    exportToCsv("Delay-Discounting-Task_" + Site + "_" + SubjectID + "_" + Session + "_" + moment().format() + ".csv",
		  			jsPsych.data.getData());
			},'json').error(function() {
			    exportToCsv("Delay-Discounting-Task_" + Site + "_" + SubjectID + "_" + Session + "_" + moment().format() + ".csv",
		  			jsPsych.data.getData());			
	    });
	     
       }
   });
			  
function makeUnique( data, prefix ) {
    
    var build, key, destKey, value;
    
    build = {};
    if (typeof data === "object") {
	if (data instanceof Array) { // don't change the array, only traverse
	    for (var i = 0; i < data.length; i++) {
		data[i] = makeUnique(data[i], prefix);
	    }
	    return data;
	} else {
	    for (key in data) {
		// Get the destination key
		destKey = prefix + key;

		// Get the value
		value = data[key];

		value = makeUnique(value, prefix);

		build[destKey] = value;
	    }
	}	
    } else {
	return data;
    }

    return build;
}
    
</script>
</html>
    
