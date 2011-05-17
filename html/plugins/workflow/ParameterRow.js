dojo.provide( "plugins.workflow.ParameterRow");

// HAS A 
dojo.require("plugins.upload.FileUpload");

// INHERITS
dojo.require("plugins.core.Common");

dojo.declare( "plugins.workflow.ParameterRow",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "workflow/templates/parametersrow.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

validInput : true,	// FILE PRESENT BOOLEAN

// CORE WORKFLOW OBJECTS
core : null,

/////}

constructor : function(args) {
	this.passedArgs = args;
	////console.dir(this.passedArgs);

	this.core = args.core;

	for ( var key in this.passedArgs )
	{
	}

	// POPULATE this.parameterObject TO BE USED IN saveStageParameter
	this.parameterObject = new Object;
	for ( var key in args )
	{
		if ( key == "core" )	continue;
		this.parameterObject[key] = args[key];
	}

},

postCreate : function() {
	this.startup();
},

startup : function () {

	this.inherited(arguments);

	// SET TOGGLE VISIBILITY WITH name ELEMENT ONCLICK
	this.setToggle();

	// SET DEFAULT INPUT VALUES (I.E., description AND value TO required|optional)
	//this.setDefaults();

	// SET LISTENERS FOR EDIT ONCLICK
	this.setEditOnclicks();

	// SET this.nameNode.parentWidget FOR 'ONCLICK REMOVE' IN Workflow.setDropTarget
	this.nameNode.parentWidget = this;

	// SET BROWSE BUTTON IF FILE OR DIRECTORY
	this.setBrowseButton();

	// SET UPLOAD BUTTON IF INPUT FILE ONLY
	if ( this.valuetype.match(/^file/) && this.paramtype == "input" )
		this.setFileUpload();

	// SET UPLOAD BUTTON ONCLICK LISTENER
	if ( this.valuetype.match(/^file/) )
		this.setFileDownload();

	if ( this.valuetype.match(/^flag$/) )
		this.setCheckbox();	
},


getInputValue : function () {
// GET INPUTS FROM THE FORM
	var inputValue;
	if ( this.valueNode.innerHTML )
	{
		inputValue = this.valueNode.innerHTML;
	}

	if ( inputValue && inputValue.match(/^<input type="checkbox">/) )
	{
		inputValue = this.valueNode.firstChild.checked;
	}

	var name = "valueNode";
	var value;

	//// WIDGET COMBO BOX
	//if ( this.valueNode.getValue != null )
	//{
	//	value = this.valueNode.getValue();
	//}
	//
	//// TEXT INPUT OR HTML COMBO BOX
	//else
	//{
	//	value = this.valueNode.value;

		if ( this.formInputs[name] == "word" )
		{
			value = this.cleanWord(value);
		}
		else if ( this.formInputs[name] == "phrase" )
		{
			value = this.cleanEnds(value);
		}

	//	this.valueNode.value = value;
	//}

	return inputValue;
},



destroy: function (preserveDom) {
// overridden destroy method. call the parent method


	this.inherited(arguments);
},


/////////////////		FLAG METHODS         
setCheckbox : function () {
// ADD A CHECKBOX TO THE valueNode IF THE INPUT TYPE IS A FLAG

	var checkbox = document.createElement('input');
	checkbox.type = "checkbox";
	if ( this.value == "checked" )
		checkbox.checked = true;
	else
		checkbox.checked = false;
	this.valueNode.innerHTML = '';
	this.valueNode.appendChild(checkbox);

	// SET ONCHANGE LISTENER
	var thisObject = this;
	dojo.connect(checkbox, "onchange", dojo.hitch(this, function(event)
		{

			//Stop Event Bubbling
			event.stopPropagation(); 

			// GET INPUTS
			var inputs = thisObject.getEditedInputs(thisObject);
			if ( inputs == null ) return;

			// SAVE STAGE PARAMETER
			thisObject.saveStageParameter(inputs);

			// UPDATE validInputs IN Parameter AND StageRow
			thisObject.core.parameters.checkValidInputs();
		}
	));
},


/////////////////		BROWSE METHODS         
setBrowseButton : function() {

	// RETURN IF THE PARAMTYPE IS NOT file* OR director*
	if (! this.valuetype.match(/^file/) && ! this.valuetype.match(/^director/))	return;

	dojo.connect(this.browseButton, "onclick", this, dojo.hitch( this, function(event)
		{
			// DO FILTER REPORT
			this.openFileManager();
		}
	));
},

openFileManager : function() {
/* OPEN FILE MANAGER TO ALLOW SELECTION OF FILE AS ARGUMENT VALUE.
	PASS THIS PARAMETER ROW OBJECT AS A PARAMETER TO BE USED IN
	THE CALLBACK TO SET THE ARGUMENT VALUE (LOCALLY AND REMOTELY).
*/

	if ( this.core.fileManager == null )
	{
		return;
	}

	if ( this.paramtype != "input" )
	{
		this.core.fileManager.disableMenus();
	}
	else {
		this.core.fileManager.enableMenus();	
	}

	////console.dir(this);

	this.core.fileManager.show(this);
},


/////////////////		CHECKINPUT METHODS         
checkEditedInput : function(node, inputValue, force) {

	this.checkInput(node, inputValue, force);

	// REFRESH PARENT INFO PANE'S this.validInputs
	// UPDATE this.validInputs IN PARENT WIDGET
	this.core.parameters.checkValidInputs();		

},

checkInput : function (node, inputValue, force) {
// CHECK IF INPUT IS VALID OR IF FILE/DIRECTORY IS PRESENT.
// SET NODE CSS ACCORDINGLY.


	// NO NEED TO CHECK FLAG
	if ( this.valuetype == "flag" )
	{
		this.validInput == true;	
		this.setValid(node);
		return;
	}

	// IF EMPTY, CHECK IF REQUIRED OR ESSENTIAL
	if ( inputValue == null || inputValue == '' )
	{
		// IF ESSENTIAL/REQUIRED, SET AS INVALID
		if ( this.discretion == "essential"
			||	this.discretion == "required" )
		{
			// SET this.validInput AS FALSE AND SET invalid CSS
			this.setInvalid(node);
			this.validInput = false;

			// SET STAGE PARAMETER'S  AS false
			Agua.setParameterValidity(this, false);

			//// SET AS REQUIRED
			//this.inputRequired(node);
		}
		else
		{
			// SET STAGE PARAMETER  AS TRUE
			Agua.setParameterValidity(this, true);

			// SET this.validInput AS TRUE
			this.setValid(node);
			this.validInput = true;

			// REMOVE ANY 'REQUIRED' CSS, E.G., FILE NOT DIRECTORY
			dojo.removeClass(node, 'required');
		}

		// SET PARENT INFO PANE'S validInputs FLAG AND
		// CALL stages.updateValidity
		this.core.parameters.checkValidInputs();
	}

	// DO TEXT INPUT IF TYPE IS NOT file, files, directory, OR directories
	else if ( ! this.valuetype.match(/^file$/)
			&& ! this.valuetype.match(/^directory$/) )
	{
		this.checkTextInput(node, inputValue);
		return;
	}

	// OTHERWISE, DO FILE INPUT
	else
	{
		force = false;
		this.checkFile(node, inputValue, force);
		return;
	}
},

checkTextInput : function (node, inputValue) {
/* CHECK IF THE VALUE IS NULL OR INVALID.
// IF SO, CHECK IF VALUE REQUIRED AND LABEL ACCORDINGLY */

	// ADD invalid CLASS IF INPUT IS NOT VALID
	this.validInput = this.isValidInput(this.invalidInputs, "valueNode", inputValue);

	// IF THE INPUT IS NOT VALID, ADD THE invalid CLASS.
	// ADD THE required CLASS IF 
	if ( this.validInput == false
		&& this.paramtype == "input" )
	{
		this.validInput = false;
		this.setInvalid(node);

		// SET STAGE PARAMETER'S  AS true
		Agua.setParameterValidity(this, false);
	}

	// IF THE INPUT IS VALID, TOGGLE ITS 'required|satisfied'
	// ACCORDINGLY IF THE INPUT IS REQUIRED. IF NOT, DO NOTHING
	else
	{
		// SET STAGE PARAMETER'S  AS true
		Agua.setParameterValidity(this, true);
		this.setValid(node);
	}


	//// MAKE THE PARENT Parameters CHECK FOR VALID INPUTS
	//// AMONG ALL OF ITS ParameterRow CHILDREN AND UPDATE
	//// THE DISPLAY OF THE STAGE IN THE WORKFLOW ACCORDINGLY.
	//// ALSO SET PARENT INFO PANE'S validInputs FLAG AND
	//// CALL stages.updateValidity
	this.core.parameters.checkValidInputs();

},	//	checkTextInput


checkFile : function (node, inputValue, force) {
/* CHECK IF FILE/DIRECTORY IS PRESENT ON SERVER.
// CALL PARENT WIDGET TO UPDATE ITS validInputs SLOT.
// NB: RETURN NULL IF inputValue IS EMPTY OR NULL */

	if ( inputValue == null || inputValue == '' )	return null;

	// GENERATE QUERY JSON FOR THIS WORKFLOW IN THIS PROJECT
	var query = new Object;

	// SET requestor = THIS_USER IF core.parameters.shared IS TRUE
	if ( this.core.parameters.shared == true )
	{
		query.username = this.username;
		query.requestor = Agua.cookie('username');
	}
	else
	{
		query.username = Agua.cookie('username');
	}

	query.sessionId = Agua.cookie('sessionId');
	query.project = this.project;
	query.workflow = this.workflow;
	query.mode = "checkFile";
	query.filepath = inputValue;

	var url = this.randomiseUrl(Agua.cgiUrl + "workflow.cgi");

	// ADD RANDOM NUMBER CONTENT TO DISAMBIGUATE xhrPut REQUESTS ON SERVER
	var content = Math.floor(Math.random()*1000000000000);

	// SEND TO SERVER
	var thisObject = this;
	dojo.xhrPut(
		{
			url: url,
			contentType: "text",
			preventCache: true,
			sync : false,
			handleAs: "json",
			content: content,
			putData: dojo.toJson(query),
			timeout: 20000,
			load: function(response, ioArgs) {
				thisObject.handleCheckfile(thisObject.valueNode, response);
			},
			error: function(response, ioArgs) {

			}
		}
	);	
},


handleCheckfile : function (node, fileinfo) {
// HANDLE JSON RESPONSE FROM checkFile QUERY


	if ( fileinfo == null || fileinfo.exists == null )
	{
		return;
	}

	//if ( this.discretion == "required" )
	//{
	//	this.validInput = true;
	//	this.inputSatisfied(node);
	//
	//	// SET STAGE PARAMETER'S isValid AS true
	//	Agua.setParameterValidity(this, true);
	//	return;
	//}
	//

	// IF THE FILE EXISTS, SET this.validInput TO TRUE
	// AND ADD filePresent AND inputSatisfied CSS CLASSES 
	if ( fileinfo.exists == "true" )
	{
		// SET FILE PRESENT CSS
		this.filePresent(node);

		// IF file, files OR directory, directories SET CSS TO satisfied
		if( (this.valuetype == "directory" && fileinfo.type == "directory")
		   || (this.valuetype == "file" && fileinfo.type == "file") )
		{

			// SET SATIFIED/REQUIRED CSS AND this.validInput			
			this.setValid(node);
			this.validInput = true;

			//// SET INPUT SATISFIED CSS
			//if ( this.discretion == "essential" || this.discretion == "required" )
			//{
			//	this.inputSatisfied(node);
			//}

			// SET STAGE PARAMETER'S isValid AS true
			Agua.setParameterValidity(this, true);
		}

		// OTHERWISE, ITS A DIRECTORY WHEN A FILE IS REQUIRED, OR VICE-VERSA
		// SO SET CSS CLASS TO required
		else
		{
			this.setInvalid(node);
			this.validInput = false;
			//this.inputRequired(node);

			// SET STAGE PARAMETER'S isValid AS false
			Agua.setParameterValidity(this, false);
		}
	}	// fileinfo.exists == true

	// OTHERWISE, ADD fileMissing AND inputRequired CSS
	// CLASSES AND SET this.validInput TO FALSE
	else
	{

		// SET FILE MISSING CSS
		this.fileMissing(node);

		// IF FILE MUST BE PHYSICALLY PRESENT (I.E., IT'S essential)
		// SET required CSS AND this.validInput TO FALSE			
		if ( this.discretion == "essential" )
		{
			this.validInput = false;
			this.setInvalid(node);
			//this.inputRequired(node);

			// SET STAGE PARAMETER'S isValid AS false
			Agua.setParameterValidity(this, false);
		}

		// FILE IS NOT REQUIRED TO BE PHYSICALLY PRESENT
		// (I.E., IT'S NOT essential).
		// SET satisfied CSS AND this.validInput TO TRUE
		else
		{
			this.validInput = true;
			//this.inputSatisfied(node);
			this.setValid(node);

			// SET STAGE PARAMETER'S isValid AS true
			Agua.setParameterValidity(this, true);
		}
	}



	// MAKE PARENT WIDGET CHECK ALL INPUTS ARE VALID AND SET
	// ITS OWN isValid FLAG AND CSS ACCORDINGLY



	// SET PARENT INFO PANE'S validInputs FLAG AND
	// CALL stages.updateValidity
	this.core.parameters.checkValidInputs();
},



isValidInput : function (invalidInputs, name, value) {
/* RETURN true IF THE NAME DOES NOT HAVE AN ENTRY IN
// invalidInputs (A STRING OR ARRAY OF INVALID VALUES).
// ALSO CHECK FOR INTEGERS */

	if ( invalidInputs[name] == null ) return true;

	// CHECK FOR NUMBER VALIDITY		
	if ( this.valuetype == "integer" )
	{
		if ( value != null
				&& value != ''
				&& ! value.match(/^\s*[\d\.]+\s*$/) )
		{
			return false;
		}
	}
	else return true;

	//// ITS OPTIONAL AND EMPTY, RETURN IS VALID
	//if ( this.discretion == "optional"
	//		&& value == '' )
	//{
	//	return true;
	//}

	//// OTHERWISE, CHECK FOR TEXT VALIDITY
	//for ( var i = 0; i < invalidInputs[name].length; i++ )
	//{
	//
	//	if ( value == invalidInputs[name][i] )
	//	{
	//		return false;
	//	}
	//}

	return true;
},


filePresent : function (node) {
	dojo.removeClass(node, 'fileMissing');
	dojo.addClass(node, 'filePresent');
},

fileMissing : function (node) {
	dojo.removeClass(node, 'filePresent');
	dojo.addClass(node, 'fileMissing');
},
setInvalid : function (node) {
	dojo.addClass(node, 'invalid');
	dojo.addClass(this.domNode, 'invalid');

	if ( this.discretion == "required" || this.discretion == "essential" )
	{
		dojo.removeClass(node, 'satisfied');
		dojo.addClass(node, 'required');
		dojo.removeClass(this.domNode, 'satisfied');
		dojo.addClass(this.domNode, 'required');
	}
},

setValid : function (node) {
	dojo.removeClass(node, 'invalid');
	dojo.removeClass(this.domNode, 'invalid');

	if ( this.discretion == "essential"
		|| this.discretion == "required" )
	{
		dojo.removeClass(node, 'required');
		dojo.addClass(node, 'satisfied');
		dojo.removeClass(this.domNode, 'required');
		dojo.addClass(this.domNode, 'satisfied');
	}
},

/////////////////		TOGGLE METHODS         
setToggle : function () {
	// CONNECT TOGGLE EVENT
	var parameterRowObject = this;
	dojo.connect( this.nameNode, "onclick", function(event) {
		event.stopPropagation();
		parameterRowObject.toggle();
	});
},


toggle : function () {

	////console.dir(this);

	// TOGGLE HIDDEN TABLE
	// TO MAKE LAST ROW TAKE UP ALL OF THE REMAINING SPACE

	if ( this["hidden"].style.display == 'table' ) this["hidden"].style.display='none';
	else this["hidden"].style.display = 'table';

	// TOGGLE HIDDEN ELEMENTS
	//var array = [ "description", "notes" ];
	var array = ["descriptionNode", "typeNode", "typeTitle"];
	for ( var i in array )
	{
		if ( this[array[i]].style.display == 'table-cell' ) this[array[i]].style.display='none';
		else this[array[i]].style.display = 'table-cell';
	}

	// DO SPECIAL TOGGLE FOR UPLOAD AND BROWSE BUTTONS
	// DEPENDING ON PARAMETER TYPE: file OR directory
	var buttons;
	var end = 0;
	if ( this.paramtype == "input" )
	{
		buttons = ["browseButton", "downloadButton", "uploadButton", "fileInputMask"];
		if ( this.valuetype == "file" || this.valuetype == "files" )	end = 4;
		if ( this.valuetype == "directory" || this.valuetype == "directories" )	end = 1;
		//if ( this.valuetype == "integer" )	end = 0;
		//if ( this.valuetype == "string" )	end = 0;
		//if ( this.valuetype == "flag" )	end = 0;

		// REMOVE UPLOAD IF this.core.parameters SHARED IS TRUE
		if ( this.core.parameters.shared == true
				&& end == 4 )
		{
			end = 2;
		}		
	}
	else if ( this.paramtype == "output" )
	{
		buttons = ["browseButton", "downloadButton"];
		if ( this.valuetype == "file" ) end = 2;
	}


	// DO TOGGLE
	for ( var i = 0; i < end; i++ )
	{

		if ( this[buttons[i]].style.display == 'table-cell' ) this[buttons[i]].style.display='none';
		else this[buttons[i]].style.display = 'table-cell';
	}
},


/////////////////		EDIT VALUE METHODS         
setEditOnclicks : function () {
// ADD 'ONCLICK' EDIT VALUE LISTENERS

	if ( this.paramtype != "input" )
	{
		return;
	}
	var thisObject = this;
	var array = ["valueNode", "descriptionNode"];
	for ( var i in array )
	{
		// IGNORE IF TYPE IS FLAG (SET FLAG ONCHANGE EARLIER IN setCheckbox)
		if ( this.valuetype != "flag" )
		{

			dojo.connect(this[array[i]], "onclick", dojo.hitch(this, function(event)
				{
					var node = event.target;
					this.editParameterRow(this, node);
					event.stopPropagation(); //Stop Event Bubbling
				}
			));
		}	
	}
},

editParameterRow : function (parameterRowWidget, node) {
// EDIT INFOPANE ROW - SAVE IF 'RETURN' PRESSED,
// EXIT WITHOUT CHANGES IF 'ESCAPE' PRESSED


	// RETURN IF ALREADY EDITING PARAMETER ROW (I.E., MULTIPLE CLICKS)
	if ( this.editingParameterRow == true ) return;
	this.editingParameterRow = true;

	// GET THE TEXT VALUE
	var text = node.innerHTML;
	if ( text == null || ! text ) text = '';

	// RETURN IF THIS IS A DOUBLE-CLICK
	if ( text.match(/^<i/) ||
		text.match(/^<br/) ||
		text.match(/^<fieldset/) ||
		text.match(/^<textarea/) )
	{
		this.editingParameterRow = false;
		return;
	}

	// REMOVE THE EXISTING VALUE
	node.innerHTML = '';

	// CREATE INPUT TEXT AREA
	var textarea = document.createElement('textarea');
	textarea.rows = 1;
	dojo.addClass(textarea, 'editParameterRow');
	node.appendChild(textarea);
	textarea.value = text;
	textarea.focus();

	// SET NEW PARAMETER LISTENER
	var thisObject = this;
	dojo.connect(textarea, "onkeypress", function(evt){

		// summary: handles keyboard events
		var key = evt.charOrCode;
		evt.stopPropagation();

		if ( key == 13 )
		{
			// GET NEW VALUE OF INPUT
			var value = textarea.value;

			// REMOVE TEXTAREA
			node.removeChild(textarea);

			// SET NEW VALUE OF INPUT
			node.innerHTML = value;

			// NOT EDITING ROW ANY MORE
			thisObject.editingParameterRow = false;

			// ****** VALIDATE VALUE ******
			// IF NEW INPUT VALUE IS INVALID, ADD 'invalid' CSS CLASS AND QUIT
			var name = node.getAttribute('class').match(/^(\S+)/)[1];

			// HACK FOR value
			if ( name == "value" )	name = "valueNode";

			// CHECK FOR INVALID INPUT FOR NON-REQUIRED FIELDS
			var isValid = thisObject.isValidInput(thisObject.invalidInputs, name, value);

			// MAKE SURE A REQUIRED FIELD IS NOT EMPTY
			var required = thisObject.discretion == "required";
			var notSatisfied = 	required == true && (value == null || value == '');
			if ( isValid == false || notSatisfied )
			{
				// ADD invalid CLASS
				if ( node != null )	dojo.addClass(node, 'invalid');

				// RESTORE ORIGINAL VALUE OF INPUT
				//node.innerHTML = text;
			}
			else
			{
				// REMOVE invalid CLASS
				if ( node != null && node.removeClass )	dojo.removeClass(node, 'invalid');
			}

			// GET INPUTS
			var inputs = thisObject.getEditedInputs(parameterRowWidget);
			if ( inputs == null ) return;

			// SAVE STAGE PARAMETER
			thisObject.saveStageParameter(inputs);

			// CHECK IF THE FILE IS THERE, IF IT'S A FILE INPUT.
			// (WE ALREADY KNOW IT'S VALID IF ITS A TEXT INPUT)
			//thisObject.checkInput(node, value);
			var force = true;
			setTimeout(function(thisObj) { thisObject.checkInput(node, value); }, 100, this);

			// NO LONGER EDITING PARAMETER ROW
			thisObject.editingParameterRow = false;

		}
		else if (key == dojo.keys.ESCAPE)
		{
			// NO LONGER EDITING PARAMETER ROW
			thisObject.editingParameterRow = false;

			// REMOVE TEXTAREA
			if ( this.valuetype != "flag" )	node.removeChild(textarea);

			// RESTORE ORIGINAL VALUE
			if ( this.valuetype != "flag" )	node.innerHTML = text;	
		}
	});

	// QUIT EDIT IF FOCUS IS LOST
	dojo.connect(parameterRowWidget, "onBlur", function(evt){

		// RETURN IF editingParameterRow FLAG IS FALSE
		if ( thisObject.editingParameterRow == false )	return;

		// RETURN IF INPUT TYPE IS FLAG OR IF NODE NOT DEFINED
		if ( thisObject.type == "flag" )	return;
		if ( node == null )	return;

		// REMOVE TEXTAREA
		var textarea = node.getElementsByTagName('textarea');
		if ( textarea == null )	return;

		if ( node.removeChild )
		{
			if ( thisObject.type != "flag"
				&& thisObject.type != "select"
				&& node.removeChild )
			{
				while ( node.firstChild )
				{
					node.removeChild(node.firstChild) 	
				}
			}
		}

		// RESTORE ORIGINAL VALUE
		node.innerHTML = text;

		// UNSET editingParameterRow FLAG
		thisObject.editingParameterRow = false;
	});
},

getEditedInputs : function (widget) {
// GET INPUTS FROM THE EDITED ITEM

	var inputs = new Object;
	for ( var name in this.formInputs )
	{

		var inputbox = this.valueNode.firstChild;


		// IF ITS A FLAG, JUST GET ITS checked VALUE
		if ( this.valuetype == "flag" && name == "valueNode" )
		{
			var checkbox = this.valueNode.firstChild;
			var checked = checkbox.checked;
			var value = checked ? "checked" : "";

			inputs["value"] = value;	
		}

		// FOR SELECT BOXES, GET value
		else if ( widget[name].value )
		{
			if ( name == "valueNode" )
			{
				inputs["value"] = String(widget[name].value.toString());
			}
			else
				inputs[name] = String(widget[name].value.toString());

			// DON'T CLEAN SELECT BOX INPUT
			//inputs[name] = this.cleanWord(inputs[name]);
		}

		// FOR COMBO BOX WIDGETS, DO getValue()
		else if ( widget[name].getValue )
		{

			if ( name == "valueNode" )
				inputs["value"] = widget[name].getValue();
			else
				inputs[name] = widget[name].getValue();
		}

		// OR GET innerHTML
		else
		{

			if ( name == "valueNode" )
				inputs["value"] = widget[name].innerHTML;
			else
				inputs[name] = widget[name].innerHTML;



			// REMOVE WHITESPACE
			if ( this.formInputs[name] == "word" )
			{
				inputs[name] = this.cleanWord(inputs[name]);
			}
			else if ( this.formInputs[name] == "phrase" )
			{
				inputs[name] = this.cleanEnds(inputs[name]);	
			}

			// INSERT THE CLEANED VALUE BACK INTO THE WIDGET
			if ( name == "valueNode" )
				widget[name].innerHTML = inputs["value"];
			else
				widget[name].innerHTML = inputs[name];
		}
	}

	return inputs;
},


saveStageParameter : function (inputs) {
// UPDATE STAGE PARAMETER IN AGUA AND ON REMOTE SERVER
	////////console.dir(inputs);
	// SET savingParameter FLAG
	if ( this.savingParameter == true )
	{
		return;
	}
	this.savingParameter = true;

	// GET PREVIOUS name AND description, ETC. EDITABLE VALUES
	////////console.dir(this.passedArgs);
	////////console.dir(this.parameterObject);

	// REPLACE EXISTING VALUES IN this AND stageParameter WITH EDITED ONES
	for ( var key in inputs )
	{
		if ( key == "params" || key == "core" )	continue;
		this.passedArgs[key] = inputs[key];
		this.parameterObject[key] = inputs[key];
	}
	delete this.parameterObject.core;
	//////console.dir(this.parameterObject);

	// SET chained TO ZERO
	this.parameterObject.chained = 0;
	//////console.dir(this.parameterObject);

	//// REPLACE DEFAULT VALUES WITH ''
	//stageParameter = this.unsetDefaultInputs(stageParameter);

	// ADD STAGE PARAMETER TO Agua.stageparameters 
	// (REPLACES ANY EXISTING RECORD WITH THE SAME appname, number AND name)
	var success = Agua.addStageParameter(this.parameterObject);

	this.savingParameter = false;

}, // Parameters.saveStageParameter


/////////////////		DOWNLOAD METHODS         
setFileDownload : function(node, name) {

	dojo.connect(this.downloadButton, "onclick", this, dojo.hitch( this, function(event)
		{

			// DO FILTER REPORT
			this.downloadFile(this.valueNode.innerHTML);
		}
	));
},

downloadFile : function (filepath) {


	var query = "?mode=downloadFile";

	// SET requestor = THIS_USER IF core.parameters.shared IS TRUE
	if ( this.core.parameters.shared == true )
	{
		query += "&username=" + this.username;
		query += "&requestor=" + Agua.cookie('username');
	}
	else
	{
		query += "&username=" + Agua.cookie('username');
	}

	query += "&sessionId=" + Agua.cookie('sessionId');
	query += "&filepath=" + filepath;

	var url = Agua.cgiUrl + "download.cgi";

	var args = {
		method: "GET",
		url: url + query,
		handleAs: "json",
		timeout: 10000,
		load: this.handleDownload
	};

	// do an IFrame request to download the csv file.
	var value = dojo.io.iframe.send(args);
},

handleDownload : function (response, ioArgs) {

	if ( response.message == "ifd.getElementsByTagName(\"textarea\")[0] is undefined" )
	{
		Agua.toastMessage("Download failed: File is not present", "error");
	}	
},





/////////////////		FILEUPLOAD METHODS         
setFileUpload : function () {
/* SET FILE UPLOAD TO GENERATE UPLOAD WINDOW WHEN CLICKED
// INPUTS:
//		1. buttonNode : LISTENS FOR ONCLICK
// 		2. inputNode: VALUE WILL CHANGE TO THE UPLOADED FILE NAME */


	// SET THE PATH AS THE WORKFLOW FOLDER		
	var path = this.project + "/" + this.workflow;

	// SET HIDDEN VALUES TO GENERATE HIDDEN NODES IN FORM
	var hiddenValues = {
		username: Agua.cookie('username'),
		sessionId : Agua.cookie('sessionId'),
		path: path
	};		

	// SET CALLBACK
	var oldPath = this.valueNode.innerHTML;

	var uploadUrl = Agua.cgiUrl + "upload.cgi";
	var callback = dojo.hitch(this, function(newValue, type)
		{
			this.changeValue(this.valueNode, oldPath, newValue, "file");
		}
	);

	// BUTTONS MUST BE DEFINED
	var inputNode = this.valueNode;
	var buttonNode = this.uploadButton;


	var fileUpload = new plugins.upload.FileUpload(
		{
			path: path,
			inputNode: inputNode,
			buttonNode: buttonNode,
			url : uploadUrl,
			callback : callback,
			hiddenValues : hiddenValues
		}
	);

	// SET SIZE OF FILE INPUT
	////////console.dir(fileUpload.fileInput);
//fileUpload.fileInput.style.width = "30px";

},


changeValue : function (node, oldValue, newValue, type) {
/* 1. DISPLAY THE NEW VALUE AND ADD IT TO Agua
// 2. CHECK THE FILE IS PRESENT
// 3. UPDATE THE STAGE DISPLAY IN THE WORKFLOW
//	TO REFLECT THE STATE OF COMPLETENESS AND 
// 	VALIDITY OF ITS INPUTS  */

	// IF SOMETHING WENT WRONG, USE THE OLD VALUE
	if ( newValue == null || newValue == '' || newValue == oldValue )
	{

		// PUT THE OLD VALUE BACK IN THE TABLE
		node.innerHTML = oldValue;
		return;
	}

	// PUT THE VALUE IN THE TABLE
	node.innerHTML = newValue;

	// SAVE THIS OPTION VALUE FOR THE WORKFLOW TO THE SERVER
	var stageParameterObject = new Object;

	for ( var key in this.passedArgs )
	{
		if ( key != "core" )
		{
			stageParameterObject[key] = this.passedArgs[key];
		}
	}
	stageParameterObject.value = newValue;

	// SET USER NAME TO COMPLETE stageparameter UNIQUE KEY
	stageParameterObject.username = Agua.cookie('username');

	// ADD STAGE PARAMETER		
	Agua.addStageParameter(stageParameterObject);

	// SET FILE PRESENT CSS
	this.filePresent(node);


	// CHECK TYPE MATCHES E.G., EXPECTED DIRECTORY AND FOUND DIRECTORY
	// IF TYPE DOES NOT MATCH, SET AS INVALID AND REQUIRED
	if ( this.valuetype.substring(0,4) != type.substring(0,4) )
	{

		//// ALTHOUGH FILE EXISTS, IF ITS THE WRONG TYPE,
		// SET THIS INPUT AS VALID
		this.validInput = false;

		// SET VALID CSS
		this.setInvalid(node);

		//// SET INPUT REQUIRED CSS IF essential OR required
		//if ( this.discretion == "essential" || this.discretion == "required" )
		//	this.inputRequired(node);

		// SET STAGE PARAMETER'S isValid AS TRUE
		Agua.setParameterValidity(this, false);

	}
	else
	{
		// SET THIS INPUT AS VALID
		this.validInput = true;

		// REMOVE RED BORDER (.infopane .invalid OR .infopane .input .required)
		this.setValid(node);

		//// SET INPUT SATISFIED CSS
		//if ( this.discretion == "essential" || this.discretion == "required" )
		//	this.inputSatisfied(node);

		// SET STAGE PARAMETER'S isValid AS TRUE
		Agua.setParameterValidity(this, true);
	}

	//////// UPDATE this.validInput AND required|satisfied NODE CSS
	////////this.checkFile(node, newValue);

	////////// UPDATE required|satisfied NODE CSS



	// UPDATE this.validInputs IN PARENT INFOPANE WIDGET
	this.core.parameters.checkValidInputs();		
},



addValue : function (node, oldValue, newValue, type) {
/* 1. DISPLAY THE NEW VALUE AND ADD IT TO Agua
// 2. CHECK THE FILE IS PRESENT
// 3. UPDATE THE STAGE DISPLAY IN THE WORKFLOW
//	TO REFLECT THE STATE OF COMPLETENESS AND 
// 	VALIDITY OF ITS INPUTS  */

	// IF SOMETHING WENT WRONG, USE THE OLD VALUE
	if ( newValue == null || newValue == '' || newValue == oldValue )
	{

		// PUT THE OLD VALUE BACK IN THE TABLE
		node.innerHTML = oldValue;
		return;
	}

	// ADD NEW VALUE TO OLD VALUE
	newValue = oldValue + "," + newValue;

	// SET INPUT VALUE TO NEW VALUE
	node.innerHTML = newValue;

	// SAVE THIS OPTION VALUE FOR THE WORKFLOW TO THE SERVER
	var stageParameterObject = new Object;
	for ( var key in this.passedArgs )
	{
		if ( key != "core" )
			stageParameterObject[key] = this.passedArgs[key];
	}
	stageParameterObject.value = newValue;

	// SET USER NAME TO COMPLETE stageparameter UNIQUE KEY
	stageParameterObject.username = Agua.cookie('username');

	// ADD STAGE PARAMETER		
	Agua.addStageParameter(stageParameterObject);

	// SET FILE PRESENT CSS
	this.filePresent(node);

	// CHECK TYPE MATCHES E.G., EXPECTED DIRECTORY AND FOUND DIRECTORY
	// IF TYPE DOES NOT MATCH, SET AS INVALID AND REQUIRED
	if ( this.valuetype != type )
		//|| ( this.valuetype == "file" && oldValue != '' )
		//|| ( this.valuetype == "directory" && oldValue != '')  )
	{

		// SET THIS INPUT AS VALID
		this.validInput = false;

		// SET VALID CSS
		this.setInvalid(node);

		//// SET INPUT REQUIRED CSS IF essential OR required
		//if ( this.discretion == "essential" || this.discretion == "required" )
		//	this.inputRequired(node);

		// SET STAGE PARAMETER'S isValid AS TRUE
		Agua.setParameterValidity(this, false);
	}
	else
	{

		// SET THIS INPUT AS VALID
		this.validInput = true;

		// REMOVE RED BORDER (.infopane .invalid OR .infopane .input .required)
		this.setValid(node);

		// SET INPUT SATISFIED CSS
		if ( this.discretion == "essential" || this.discretion == "required" )
			this.inputSatisfied(node);

		// SET STAGE PARAMETER'S isValid AS TRUE
		Agua.setParameterValidity(this, true);
	}

	//////// UPDATE this.validInput AND required|satisfied NODE CSS
	////////this.checkFile(node, newValue);

	////////// UPDATE required|satisfied NODE CSS


	// UPDATE this.validInputs IN PARENT INFOPANE WIDGET
	this.core.parameters.checkValidInputs();		
},

formInputs : {		// FORM INPUTS AND TYPES (word|phrase)
	valueNode: "word",
	descriptionNode: "phrase"
},

defaultInputs : {	// DEFAULT INPUTS
	valueNode: "",
	descriptionNode: ["Description"]
},


requiredInputs : {	// REQUIRED INPUTS CANNOT BE ''
	valueNode : 1
},

invalidInputs : {	// THESE INPUTS ARE INVALID

	valueNode : ["essential", "required", "optional"]
	//description: "Description",
	//notes: "Notes"
}



});
