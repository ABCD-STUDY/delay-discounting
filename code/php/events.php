<?php

  session_start(); /// initialize session, we will use session variables to store the subjid

  include($_SERVER["DOCUMENT_ROOT"]."/code/php/AC.php");
  $user_name = check_logged(); /// function checks if user is logged in

  if (!$user_name || $user_name == "") {
     echo (json_encode ( array( "message" => "no user name" ) ) );
     return; // nothing
  }

  $permissions = list_permissions_for_user( $user_name );

  // find the first permission that corresponds to a site
  // Assumption here is that a user can only add assessment for the first site he has permissions for!
  $site = "";
  foreach ($permissions as $per) {
     $a = explode("Site", $per); // permissions should be structured as "Site<site name>"

     if (count($a) > 0) {
        $site = $a[1];
	break;
     }
  }
  if ($site == "") {
     echo (json_encode ( array( "ok" => 0, "message" => "Error: no site assigned to this user" ) ) );
     return;
  }

  // Both the subject id and the visit (session) are used to make the assessment unique
   $subjid = "";
   $sessionid = "";
   $active_substances = array();
   $run = "";
   if ( isset($_SESSION['ABCD']) && isset($_SESSION['ABCD']['delay-discounting']) ) {
      if (isset($_SESSION['ABCD']['delay-discounting']['subjid'])) {  
         $subjid  = $_SESSION['ABCD']['delay-discounting']['subjid'];
      }
      if (isset($_SESSION['ABCD']['delay-discounting']['sessionid'])) {
         $sessionid  = $_SESSION['ABCD']['delay-discounting']['sessionid'];
      }
      if (isset($_SESSION['ABCD']['delay-discounting']['run'])) {
         $run  = $_SESSION['ABCD']['delay-discounting']['run'];
      }      
   }
   if ($subjid == "") {
     echo(json_encode ( array( "ok" => 0, "message" => "Error: no subject id assigned" ) ) );
     return;
   }
   if ($sessionid == "") {
     echo(json_encode ( array( "ok" => 0, "message" => "Error: no session specified" ) ) );
     return;
   }
   if ($run == "") {
     echo(json_encode ( array( "ok" => 0, "message" => "Error: no run specified" ) ) );
     return;
   }
  $action = "save";
  if (isset($_POST['action'])) {
    $action = $_POST['action'];
  }

  // this event will be saved at this location
  $events_file = $_SERVER['DOCUMENT_ROOT']."/applications/delay-discount/data/" . $site . "/ded_".$subjid."_".$sessionid."_".$run.".json";

  if ($action == "test") {
     // test if the current file exists already
     if (file_exists($events_file)) {
       echo(json_encode ( array( "message" => "Error: this session already exists, overwrite session is not possible", "ok" => "0" ) ) );
       return;
     }
     echo(json_encode ( array( "message" => "save ok", "ok" => "1" ) ) );
     return;
  } else if ($action == "mark" ){
      // for now ignore this
      return;
  } else {
      // we make the assumption this is now action == "save"
  }

  $dd = $_SERVER['DOCUMENT_ROOT']."/applications/delay-discount/data/" . $site;
  if (!file_exists($dd)) {
     mkdir($dd,0777);
  }

  if (file_exists($events_file)) {
     echo(json_encode ( array( "ok" => 0, "message" => "Error: this session already exists, overwrite session is not possible" ) ) );
     return;
  }
  
  $ar = array( "data"              => [],
               "ded_server_date"    => date("Y/m/d"),
               "ded_server_time"    => date("h:i:sa"),
	       "ded_site"          => $site,
	       "ded_subject_id"    => $subjid,
	       "ded_run"           => $run,
	       "ded_event_name"    => $sessionid,
	       "id_redcap"         => $subjid,
	       "redcap_event_name" => $sessionid
  );
  if (isset($_POST['toplevel'])) {
       foreach($_POST['toplevel'] as $key => $value) {
   	  $ar[$key] = $value;
       }
  }
  if (isset($_POST['data'])) {
     $ar['data'] = json_decode($_POST['data'], true);
  }
  if (isset($_POST['date'])) {
     $ar['ded_assessment_date'] = $_POST['date'];
  }
  file_put_contents($events_file, json_encode( $ar, JSON_PRETTY_PRINT ));
  echo (json_encode( array( "ok" => "1", "message" => "Saved session" ) ));
?>
