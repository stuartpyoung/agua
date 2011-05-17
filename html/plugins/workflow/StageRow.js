dojo.provide( "plugins.workflow.StageRow");

dojo.require("plugins.core.Common");

dojo.declare( "plugins.workflow.StageRow",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	/////}

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "workflow/templates/stagerow.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// OR USE @import IN HTML TEMPLATE
cssFiles : [ dojo.moduleUrl("plugins") + "/workflow/css/stagerow.css" ],

// PARENT plugins.workflow.Apps WIDGET
parentWidget : null,

// APPLICATION OBJECT
application : null,

// CORE WORKFLOW OBJECTS
core : null,

constructor : function(args) {

	this.loadCSS();

	this.core = args.core;

	this.parentWidget = args.parentWidget;

	this.application = new Object;
	for ( var key in args )
	{
		if ( key != "parentWidget" )
		{
			this.application[key] = args[key];
		}
	}

	//this.inherited(arguments);
},

// RETURN A COPY OF this.application
getApplication : function () {
	return dojo.clone(this.application);
},

// SET this.application TO THE SUPPLIED APPLICATION OBJECT
setApplication : function (application) {
	this.application = application;

	return this.application;
},

// SET THE NUMBER NODE TO THE stage.number 
setNumber : function (number) {

	this.application.number = number;
	this.application.appnumber = number;
	this.numberNode.innerHTML = number;
},


postCreate : function() {

	this.startup();
},

/* SET this.name.parentWidget = this

	FOR RETRIEVAL OF this.application WHEN MENU IS CLICKED

	REM: remove ONCLICK BUBBLES ON stageRow.name NODE RATHER THAN ON node. 

	I.E., CONTRARY TO DESIRED, this.name IS THE TARGET INSTEAD OF THE node.

	ALSO ADDED node.parentWidget = stageRow IN Stages.setDropTarget()

*/
startup : function () {

	this.inherited(arguments);

	this.name.parentWidget = this;

	this.setNumber(this.application.number);
},


/*	CHECK ALL THE PARAMETERS HAVE VALID INPUTS AND CHANGE CSS ACCORDINGLY.

 	RETURN VALUE OF validInputs AS true IF ALL PARAMETERS ARE SATISFIED.

 	OTHERWISE, RETURN false.

	PROCESS FOR CHECKING VALIDITY FOR EACH PARAMETER

		1. SET this.validInputs TO TRUE, UPDATE ALONG THE WAY

		2. CHECK ONLY inputs (IGNORE outputs AND resources).

			FOR EACH INPUT:

			2.1 CHECK IF VALIDITY HAS ALREADY BEEN COMPUTED AND STORED

				IN Agua.getParameterValidity BOOLEAN. USE IF AVAILABLE AND NEXT 

			2.2 OTHERWISE, SET Agua.setParameterValidity FOR EACH

				PARAMETER AS FOLLOWS:

				a. IF FILE/DIR, IGNORE IF OPTIONAL UNLESS FILEPATH SPECIFIED

				b. IF NON-OPTIONAL, ADD TO ARRAY OF files TO CHECK

				c. IGNORE FLAGS, CHECK INTS AND NON-OPTIONAL TEXT

			2.3 FOR 'UNKNOWN' FILES/DIRS, DO A BATCH FILE CHECK ON

				THE REMOTE SERVER. SET validInputs TO FALSE IF ANY

				FILES/DIRS ARE MISSING

		3. ADJUST CSS CLASS OF StageRow ACCORDING TO VALUE OF validInputs

		4. RETURN BOOLEAN VALUE OF validInputs

*/

checkValidParameters : function (force) {
	//////console.dir(this);

	// SET this.isValid TO DEFAULT true
	this.isValid = true;

	// GET STAGE PARAMETERS
	var stageParameters = Agua.getStageParameters(this.application);
	if ( stageParameters == null )
	{
		return false;
	}

	// GET ALL REQUIRED/ESSENTIAL INPUT FILE/DIRECTORY PARAMETERS
	var files = new Array;
	for ( var i = 0; i < stageParameters.length; i++ )
	{
		if ( stageParameters[i].paramtype != "input" )
		{
			continue;
		}

		stageParameters[i].value = stageParameters[i].value.replace(/%project%/, stageParameters[i].project);
		stageParameters[i].value = stageParameters[i].value.replace(/%workflow%/, stageParameters[i].workflow);


		// UNLESS force SPECIFIED, GET VALIDITY IF EXISTS
		if ( force != null )
		{
			var isValid = Agua.getParameterValidity(stageParameters[i]);

			// IF isValid EXISTS AND IS true, MOVE TO NEXT STAGE PARAMETER
			if ( isValid != null && ! force )
			{
				// SET this.isValid TO FALSE IF THIS STAGE PARAMETER isValid IS FALSE
				if (isValid == false )	this.isValid = false;
				continue;
			}
		}

		// SAVE UNKNOWN FILE/DIRECTORY FOR checkFiles LATER ON
		if ( stageParameters[i].valuetype.match(/^file$/)
			||  stageParameters[i].valuetype.match(/^directory$/) )
		{
			var filepath = stageParameters[i].value;

			// CHECK NON-OPTIONAL FILEPATHS
			// required: FILE PATH MUST BE NON-EMPTY
			// essential: FILE/DIRECTORY MUST BE PHYSICALLY PRESENT
			if ( stageParameters[i].discretion == "required"
				|| stageParameters[i].discretion == "essential" )
			{
				var exists = stageParameters[i].exists;

				// SET this.isValid = false IF FILE/DIR IS NULL OR EMPTY
				if ( filepath == null || filepath == '' )
				{
					this.isValid = false;
					Agua.setParameterValidity(stageParameters[i], false);
				}

				// OTHERWISE, IF FILE/DIR IS SPECIFIED BUT NOT KNOWN TO EXIST,
				// PUSH ONTO files ARRAY
				// IGNORE IF CHAINED IS DEFINED AND 1
				else if ( (stageParameters[i].exists == null
					|| stageParameters[i].exists == false)
					&& ( stageParameters[i].chained == null
					|| ! stageParameters[i].chained ))
				{
					files.push(stageParameters[i]);
				}
			}

			// THIS IS AN OPTIONAL PARAMETER SO ITS VALID IF EMPTY BUT
			// IF ITS NOT EMPTY, CHECK IT EXISTS
			else 
			{
				// IF EMPTY 
				if ( stageParameters[i].value == null
					|| stageParameters[i].value == '' )
				{
					Agua.setParameterValidity(stageParameters[i], true);
				}
				else
				{
					Agua.setParameterValidity(stageParameters[i], true);
					files.push(stageParameters[i]);
				}
			}
		}

		// FLAGS ARE AUTOMATICALLY VALID
		else if ( stageParameters[i].valuetype.match(/flag/) )
		{
			Agua.setParameterValidity(stageParameters[i], true);
		}

		// CHECK INTEGERS
		else if ( stageParameters[i].valuetype.match(/integer/) )
		{
			if ( stageParameters[i].discretion != "optional" )
			{
				// SET EMPTY NON-OPTIONAL INTEGER AS FALSE
				if ( stageParameters[i].value == null
						|| stageParameters[i].value == '' )
				{
					this.isValid = false;
					Agua.setParameterValidity(stageParameters[i], false);
				}
				// NON-OPTIONAL INTEGER BUT IS NOT A PROPER NUMBER
				// ////SO SET TO false
				else if (! stageParameters[i].value.match(/^\s*[\d\.]+\s*$/) )
				{
					this.isValid = false;
					Agua.setParameterValidity(stageParameters[i], false);
				}
				// OTHERWISE, ITS A CORRECT OPTIONAL INTEGER SO SET TO true
				else
				{
					Agua.setParameterValidity(stageParameters[i], true);
				}
			}
			else
			{
				// SET OPTIONAL INTEGER TO false IF ITS NON-EMPTY BUT NOT
				// AN INTEGER
				if ( stageParameters[i].value != null
						&& stageParameters[i].value != ''
						&& ! stageParameters[i].value.match(/^\s*[\d\.]+\s*$/) )
				{
					this.isValid = false;
					Agua.setParameterValidity(stageParameters[i], false);
				}

				// OTHERWISE, ITS EITHER EMPTY OR AN INTEGER SO SET TO true
				else
				{
					Agua.setParameterValidity(stageParameters[i], true);
				}
			}
		}

		// CHECK TEXT INPUTS
		else
		{
			if ( stageParameters[i].discretion != "optional" )
			{
				if ( stageParameters[i].value == null
					|| stageParameters[i].value == '' )
				{
					this.isValid = false;
					Agua.setParameterValidity(stageParameters[i], false);
				}
				else
				{
					Agua.setParameterValidity(stageParameters[i], true);
				}
			}
			else
			{
				// THIS IS AN OPTIONAL PARAMETER SO, EMPTY OR NOT, ITS VALID
				Agua.setParameterValidity(stageParameters[i], true);
			}
		}
	}

	if ( files.length == 0 )
	{

		if ( this.isValid == false || this.isValid == null )
			this.setInvalid();
		else this.setValid();

		return this.isValid;
	}



	// GET FILEINFO FROM REMOTE FILE SYSTEM
	var url = Agua.cgiUrl + "workflow.cgi";
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.project = this.application.project;
	query.workflow = this.application.workflow;
	query.mode = "checkFiles";
	query.files = files;

	// SEND TO SERVER
	var thisObject = this;
	var xhrputReturn = dojo.xhrPut(
		{
			url: url,
			contentType: "text",
			sync : true,
			handleAs: "json",
			putData: dojo.toJson(query),
			timeout: 20000,
			load: function(fileinfos, ioArgs) {
				if ( fileinfos == null
					|| fileinfos.length == null
					|| ! fileinfos.length )	return;

				for ( var i = 0; i < fileinfos.length; i++ )
				{
					if ( fileinfos[i].exists == "true" )
					{
						// DO Agua.setFileExists
						//var success = Agua.setFileExists(files[i], true);

						// DO Agua.setParameterValidity
						Agua.setParameterValidity(files[i], true);
					}
					else
					{
						if ( files[i].discretion == "essential" )
						{


							thisObject.isValid = false;

							// DO Agua.setParameterValidity
							Agua.setParameterValidity(files[i], false);
						}
						else
						{
							// DO Agua.setParameterValidity
							Agua.setParameterValidity(files[i], true);
						}
					}
				}

			},
			error: function(response, ioArgs) {

			}
		}
	);	

	//////console.dir(xhrputReturn);

	//var fileinfos = this.getPutResult(url, query);
	//
	//if ( fileinfos == null
	//	|| fileinfos.length == null
	//	|| ! fileinfos.length )	return;
	//	
	//for ( var i = 0; i < fileinfos.length; i++ )
	//{
	//	if ( fileinfos[i].exists == "true" )
	//	{
	//		// DO Agua.setFileExists
	//		var success = Agua.setFileExists(files[i], true);
	//
	//		// DO Agua.setParameterValidity
	//		Agua.setParameterValidity(files[i], true);
	//
	//	}
	//	else
	//	{
	//		if ( files[i].discretion == "essential" )
	//		{
	//			this.isValid = false;
	//
	//			// DO Agua.setParameterValidity
	//			Agua.setParameterValidity(files[i], false);
	//		}
	//		else
	//		{
	//			// DO Agua.setParameterValidity
	//			Agua.setParameterValidity(files[i], true);
	//		}
	//	}
	//}

	if ( this.isValid == false || this.isValid == null ) this.setInvalid();
	else this.setValid();

	return this.isValid;
},


toggle : function () {

	var array = [ "executor", "localonly", "location", "description", "notes" ];
	for ( var i in array )
	{
		if ( this[array[i]].style.display == 'table-cell' ) this[array[i]].style.display='none';
		else this[array[i]].style.display = 'table-cell';
	}
},

setValid : function () {
	//for ( var key in this.core )
	//{
	//}

	dojo.removeClass(this.domNode, 'unsatisfied');
	dojo.addClass(this.domNode, 'satisfied');

	this.isValid = true;
	var stagesWidget = this.core.stages;
	stagesWidget.updateValidity();	
},

setInvalid : function () {
	for ( var key in this.core )
	{
	}

	dojo.removeClass(this.domNode, 'satisfied');
	dojo.addClass(this.domNode, 'unsatisfied');

	this.isValid = false;
	var stagesWidget = this.core.stages;
	stagesWidget.updateValidity();	
}


/* SAVE THIS FOR FURTHER TESTING:
getPutResult : function (url, query) {

	var xreturn = dojo.xhrPut(
		{
			url: url,
			contentType: "text",
			sync : true,
			handleAs: "json",
			putData: dojo.toJson(query),
			timeout: 20000,
			load: function(fileinfos, ioArgs) {
			},
			error: function(response, ioArgs) {
			}
		}
	);	
	////console.dir(xreturn);

	var putResult = dojo._contentHandlers[xreturn.ioArgs.handleAs](xreturn.ioArgs.xhr);

	return putResult;
	//return xreturn.ioArgs.xhr.responseText;
}
*/



});
