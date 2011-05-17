dojo.provide("plugins.admin.Parameters");

// ALLOW THE USER TO ADD, REMOVE AND MODIFY PARAMETERS

//dojo.require("dijit.dijit"); // optimize: load dijit layer
dojo.require("dijit.form.Button");
dojo.require("dijit.form.TextBox");
dojo.require("dijit.form.Textarea");
dojo.require("dojo.parser");
dojo.require("dojo.dnd.Source");
dojo.require("plugins.core.Common");

// HAS A
dojo.require("plugins.admin.ParameterRow");

dojo.declare(
    "plugins.admin.Parameters",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	////}

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "admin/templates/parameters.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

//addingParameter STATE
addingParameter : false,

// OR USE @import IN HTML TEMPLATE
cssFiles : [ "plugins/admin/css/parameters.css" ],

// PARENT WIDGET
parentWidget : null,

constructor : function(args) {
	// GET INFO FROM ARGS
	this.parentWidget = args.parentWidget;

	// LOAD SORIA AND FILEPICKER CSS
	this.loadCSS();		
},

postCreate : function() {

	this.startup();
},


startup : function () {

	// COMPLETE CONSTRUCTION OF OBJECT
	this.inherited(arguments);	 

	// ADD TO TAB CONTAINER		
	this.tabContainer.addChild(this.parametersTab);
	this.tabContainer.selectChild(this.parametersTab);

	// SET LISTENER FOR PARAM TYPES COMBO
	this.setParamTypesCombo();

	// SET PARAMETERS COMBO
	this.setAppsCombo();

	// SET NEW PARAMETER FORM
	this.setForm();

	// SET TRASH DROP TARGET
	this.setTrash();
},



reload : function () {
// RELOAD THE GROUP COMBO AND DRAG SOURCE AFTER CHANGES
// TO SOURCES OR GROUPS DATA IN OTHER TABS
	// SET PARAMTYPES COMBO
	this.setAppNamesCombo();

	// SET APPS COMBO
	this.setAppsCombo();

	// SET DRAG SOURCE
	this.setDragSource();
},

setAppsCombo : function (type) {
// SET PARAMETERS COMBO BOX

	// GET PARAMETERS NAMES		
	var apps = dojo.clone(Agua.apps);

	var typesHash = new Object;
	for ( var i in apps )
	{
		if ( typesHash[apps[i]] == null )
		{
			typesHash[apps[i].type] = 1;
		}
	}
	var typesArray = this.hashkeysToArray(typesHash);

	typesArray = this.sortNoCase(typesArray);
	typesArray.splice(0,0, 'Order by Type');
	typesArray.splice(0,0, 'Order by Name');


	// SET STORE
	var data = {identifier: "name", items: []};
	for ( var i = 0; i < typesArray.length; i++ )
	{
		data.items[i] = { name: typesArray[i]	};
	}
	var store = new dojo.data.ItemFileWriteStore(	{	data: data	}	);

	this.appsCombo.popupClass = "parameters appsCombo dijitReset dijitMenu";
	this.appsCombo.wrapperClass = "parameters dijitPopup";
	this.appsCombo.itemHeight = 30;


	// SET COMBO
	this.appsCombo.store = store;
	this.appsCombo.startup();

	// SET COMBO VALUE
	var firstValue = typesArray[0];
	this.appsCombo.setValue(firstValue);

	// CONNECT ONCLICK WITH dojo.connect TO BUILD TABLE
	var appsObject = this;
	dojo.connect(this.appsCombo, "onChange", function(event) {
		appsObject.setAppNamesCombo();
	});

	// SET PARAMETER NAMES COMBO
	this.setAppNamesCombo();
},


setAppNamesCombo : function () {
// SET PARAMETER NAMES COMBO DEPENDENT ON THE CURRENT SELECTION
// IN THE PARAMETERS COMBO

	// GET SOURCE ARRAY AND FILTER BY PARAMETER NAME
	var type = this.appsCombo.getValue();
	var sourceArray = dojo.clone(Agua.apps);

	// FILTER APP NAMES ARRAY BY APPS COMBO VALUE
	//sourceArray = this.filterSourceArray(sourceArray, type);
	var keyArray = ["type"];
	var valueArray = [type];
	if ( type == "Order by Name" )
		sourceArray = this.sortHasharray(sourceArray, 'name');
	else if ( type == "Order by Type" )
		sourceArray = this.sortHasharray(sourceArray, 'type');
	else
		sourceArray = this.filterByKeyValues(sourceArray, keyArray, valueArray);





	// CHECK sourceArray IS NOT NULL OR EMPTY
	if ( sourceArray == null || sourceArray.length == 0 )
	{
		return;
	}

	// SET STORE
	var data = {identifier: "name", items: []};
	for ( var i = 0; i < sourceArray.length; i++ )
	{
		data.items[i] = { name: sourceArray[i].name	};
	}
	var store = new dojo.data.ItemFileWriteStore(	{	data: data	}	);

	this.appNamesCombo.popupClass = "parameters appNamesCombo dijitReset dijitMenu";
	this.appNamesCombo.wrapperClass = "parameters dijitPopup";
	this.appNamesCombo.itemHeight = 30;


	// SET COMBO
	this.appNamesCombo.store = store;
	this.appNamesCombo.startup();

	// SET COMBO VALUE
	var firstValue = sourceArray[0].name;
	this.appNamesCombo.setValue(firstValue);

	// CONNECT ONCLICK WITH dojo.connect TO BUILD TABLE
	var appsObject = this;
	dojo.connect(this.appNamesCombo, "onChange", function(event) {
		appsObject.setDragSource();
	});

	// SET PARAMETERS COMBO
	this.setDragSource();
},




setParamTypesCombo : function () {

	// CONNECT ONCLICK WITH dojo.connect TO BUILD TABLE
	var appsObject = this;
	dojo.connect(this.paramTypesCombo, "onchange", function(event) {

		appsObject.setDragSource();
	});
},





setForm : function () {
// SET LISTENERS TO ACTIVATED SAVE BUTTON AND TO CLEAR DEFAULT TEXT
// WHEN INPUTS ARE CLICKED ON

	// SET ADD PARAMETER ONCLICK
	dojo.connect(this.addParameterButton, "onclick", dojo.hitch(this, "saveParameter", null));

	// SET ONCLICK TO CANCEL DEFAULT TEXT
	for ( var name in this.formInputs )
	{
		dojo.connect(this[name], "onclick", dojo.hitch(this, "clearValue", this[name], this.firstLetterUpperCase(name)));
		dojo.connect(this[name], "onfocus", dojo.hitch(this, "clearValue", this[name], this.firstLetterUpperCase(name)));
	}

	//// SET NEW PARAMETER LISTENER
	//var parametersObject = this;
	//dojo.connect(this.newName, "onkeypress", function(evt){
	//	var key = evt.charOrCode;
	//	if ( key == 13 )
	//	{				
	//		// SHIFT FOCUS TO NEXT INPUT
	//		parametersObject.newDescription.focus();
	//	}
	//	
	//	// REMOVE ANY LINE RETURNS
	//	this.value = this.value.replace(/\n/, '');
	//});
},
setDragSource : function () {
// SET THE DRAG SOURCE WITH PARAMETER OBJECTS

	// DELETE EXISTING TABLE CONTENT
	while ( this.dragSource.firstChild )
	{
		this.dragSource.removeChild(this.dragSource.firstChild);
	}

	// FILTER SOURCE ARRAY BY type
	var appName = this.appNamesCombo.getValue();

	var sourceArray = Agua.getParametersByAppname(appName);

	// FILTER APP NAMES ARRAY BY APPS COMBO VALUE
	var paramType = this.paramTypesCombo.value;
	var keyArray = ["paramtype"];
	var valueArray = [paramType];
	if ( paramType == "Order by Name" )
		sourceArray = this.sortHasharray(sourceArray, 'name');
	else if ( paramType == "Order by Type" )
		sourceArray = this.sortHasharray(sourceArray, 'type');
	else
		sourceArray = this.filterByKeyValues(sourceArray, keyArray, valueArray);

	sourceArray = this.sortHasharray(sourceArray, 'name');

	// CHECK IF sourceArray IS NULL
	if ( sourceArray == null )
	{
		return;
	}

	// GENERATE dataArray TO INSERT INTO DND PARAMETER TABLE
	var dataArray = new Array;
	for ( var i in sourceArray )
	{
		var data = sourceArray[i];				
		dataArray.push( { data: data, type: ["draggableItem"] } );
	}

	// GENERATE DND PARAMETER
	var dragSource = new dojo.dnd.Source(
		this.dragSource,
		{
			copyOnly: true,
			selfAccept: false,
			accept : [ "none" ]
		}
	);
	dragSource.insertNodes(false, dataArray);

	// SET TABLE ROW STYLE IN dojDndItems
	var allNodes = dragSource.getAllNodes();
	for ( var k = 0; k < allNodes.length; k++ )
	{

		// ADD CLASS FROM type TO NODE
		var node = allNodes[k];

		// SET NODE ATTRIBUTES TO IDENTIFY IT ON DND DROP AND
		// SET SOURCE HASH USED TO INSTANTIATE plugins.admin.ParameterRow
		var valueArray = [ "appname", "name", "type", "category", "argument", "value", "discretion", "format", "description", "paramtype", "args", "params", "paramFunction" ];
		var source = new Object;
		for ( var i in valueArray )
		{
			node[valueArray[i]] = source[valueArray[i]] = dataArray[k].data[valueArray[i]] ? dataArray[k].data[valueArray[i]] : this.firstLetterUpperCase(valueArray[i]);
		}
		source.parentWidget = this;

		// INSTANTIATE plugins.admin.ParameterRow AND PARAMETEREND TO NODE
		var parameterRow = new plugins.admin.ParameterRow(source);
		node.innerHTML = '';
		node.appendChild(parameterRow.domNode);
		node.widget = parameterRow;
	}

	var sourceObject = this;
	dragSource.creator = function (item, hint)
	{

		var node = dojo.doc.createElement("div");
		node.name = item.name;
		node.id = dojo.dnd.getUniqueId();
		node.className = "dojoDndItem";


		// SET FANCY FORMAT IN NODE INNERHTML
		node.innerHTML = "<table> <tr><td><strong style='color: darkred'>" + item.name + "</strong></td></tr><tr><td> " + item.description + "</td></tr></table>";

		return {node: node, data: item, type: ["text"]};
	};
},


setTrash : function () {

	var trash = new dojo.dnd.Source(
		this.trashContainer,
		{
			accept : [ "draggableItem" ]
		}
	);

	// REMOVE DUPLICATE NODES
	var thisObject = this;
	dojo.connect(trash, "onDndDrop", function(source, nodes, copy, target){
		// NODE DROPPED ON SELF --> DELETE THE NODE
		if ( target == this )
		{

			for ( var i = 0; i < nodes.length; i++ )
			{
				var node = nodes[i];

				// HACK TO AVOID THIS ERROR: node.parentNode is null
				try {
					node.parentNode.removeChild(node);


					var sourceObject = new Object;
					sourceObject.name = node.name;
					sourceObject.type = node.type;
					sourceObject.appname = node.appname;

					thisObject.deleteParameter(sourceObject);
				}
				catch (e) {
				}
			}

			// DELETE EXISTING TABLE CONTENT
			//console.dir(thisObject.trashContainer);
			while ( thisObject.trashContainer.childNodes.length > 2 )
			{

				thisObject.trashContainer.removeChild(thisObject.trashContainer.childNodes[2]);
			}


		}
	});
}, // setTrash

deleteParameter : function (sourceObject) {
// DELETE PARAMETER FROM Agua.parameters OBJECT AND IN REMOTE DATABASE

	if ( sourceObject.name == null ) 	return;
	if ( sourceObject.appname == null ) 	return;

	// REMOVING PARAMETER FROM Agua.parameters
	Agua.removeParameter(sourceObject)

	// RESET THE PARAMETERS TABLE
	this.setDragSource();

	var url = Agua.cgiUrl + "/agua?";

	// CREATE JSON QUERY
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "deleteParameter";
	query.data = sourceObject;

	// SEND TO SERVER
	dojo.xhrPut(
		{
			url: url,
			contentType: "text",
			putData: dojo.toJson(query),
			timeout: 15000,
			load: function(response, ioArgs) {
				return response;
			},
			error: function(response, ioArgs) {
				return response;
			}
		}
	);

}, // Parameters.deleteParameter


addParameter : function (sourceObject, reload) {
//	ADD PARAMETER TO Agua.parameters AND TO REMOTE DATABASE

	if ( Agua.isParameter(sourceObject.name) )
	{
		dojo.removeClass(this.newName, 'invalid');
		dojo.addClass(this.newName, 'invalid');
		this.newName.focus();
		return;
	}

	// ADD PARAMETER TO Agua.parameters ARRAY
	Agua.addParameter(sourceObject);

	// RESET PARAMETER TABLE
	this.setDragSource();

	// SAVE PARAMETER TO REMOTE DATABASE
	this.saveParameter(sourceObject, reload);

},	// Parameters.addParameter


editParameterRow : function (parameterRowWidget, node) {
// EDIT PARAMETER ROW - SAVE IF 'RETURN' PRESSED,
// EXIT WITHOUT CHANGES IF 'ESCAPE' PRESSED


	// RETURN IF ALREADY EDITING PARAMETER ROW (I.E., MULTIPLE CLICKS)
	if ( this.editingParameterRow == true ) return;
	this.editingParameterRow = true;

	// GET THE TEXT VALUE
	var text = node.innerHTML;

	node.innerHTML = '';
	if ( text == null || ! text ) text = '';

	// RETURN IF THIS IS A DOUBLE-CLICK
	if ( text.match(/^<i/) ||
		text.match(/^<br/) ||
		text.match(/^<fieldset/) ||
		text.match(/^<textarea/) )
	{
		this.editingParamterRow = false;
		return;
	}

	// CREATE INPUT TEXT AREA
	var textarea = document.createElement('textarea');
	dojo.addClass(textarea, 'editParameterRow');
	node.appendChild(textarea);
	textarea.value = text;
	textarea.focus();

	// SET NEW PARAMETER LISTENER
	var appsObject = this;
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
			appsObject.editingParameterRow = false;

			// IF NEW INPUT VALUE IS INVALID, ADD 'invalid' CSS CLASS AND QUIT
			var key = node.getAttribute('class').match(/^(\S+)/)[1];
			if ( (appsObject.isValidInput(key, value) == false
					&& appsObject.requiredInputs[key] != null) )
				//|| (appsObject.requiredInputs[key] != null
				//	&& value == null || value == '') )
			{
				dojo.addClass(node, 'invalid');

				// RESTORE ORIGINAL VALUE OF INPUT
				node.innerHTML = text;
			}

			// OTHERWISE, SAVE THE PARAMETER AND RELOAD THE DRAG SOURCE
			else
			{
				dojo.removeClass(node, 'invalid');

				// GET INPUTS
				var inputs = appsObject.getEditedInputs(parameterRowWidget);
				if ( inputs == null ) return;

				// SAVE PARAMETER
				appsObject.saveParameter(inputs, "reload");
			}
		}
		else if (key == dojo.keys.ESCAPE)
		{
			appsObject.editingParameterRow = false;

			// REMOVE TEXTAREA
			node.removeChild(textarea);

			// RESTORE ORIGINAL VALUE
			node.innerHTML = text;
		}
	});

	// QUIT EDIT IF FOCUS IS LOST
	dojo.connect(parameterRowWidget, "onBlur", function(evt){
		appsObject.editingParameterRow = false;

		// REMOVE TEXTAREA
		node.removeChild(textarea);

		// RESTORE ORIGINAL VALUE
		node.innerHTML = text;
	});

	//dojo.connect(textarea, "onblur", function(evt){
	//	appsObject.editingParameterRow = false;
	//
	//	// REMOVE TEXTAREA
	//	node.removeChild(textarea);
	//	
	//	// RESTORE ORIGINAL VALUE
	//	node.innerHTML = text;		
	//});
},
saveParameter : function (inputs, reload) {
//	SAVE A PARAMETER TO Agua.parameters AND TO REMOTE DATABASE


	if ( this.savingParameter == true )	return;
	this.savingParameter = true;

	if ( inputs == null )
	{
		inputs = this.getFormInputs(this);

		// RETURN IF INPUTS ARE NULL OR INVALID
		if ( inputs == null || this.validateFormInputs(inputs, this) == false )
		{
			this.savingParameter = false;
			return;
		}
	}

	// SET inputs APPLICATION NAME
	var appName = this.appNamesCombo.getValue();
	var appType = Agua.getAppType(appName);
	inputs.appname = appName;
	inputs.apptype = appType;

	// OTHERWISE, REMOVE ORIGINAL PARAMETER OBJECT FROM Agua.parameters 
	// THEN ADD NEW PARAMETER OBJECT TO Agua.parameters ARRAY
	//Agua.removeParameter(inputs); // not necessary as addParameter does delete first
	Agua.addParameter(inputs);

	// REDO PARAMETER TABLE
	if ( reload != null )	this.setDragSource();

	// REPLACE DEFAULTS WITH ''
	// REPLACE ' WITH "
	for ( var name in this.defaultInputs )
	{
		if ( inputs[name] != null && inputs[name] == this.defaultInputs[name] )
		{
			inputs[name] = '';
		}

		inputs[name] = inputs[name].replace(/'/g, '"');
	}

	// SAVE NEW PARAMETER TO REMOTE DATABASE
	var url = Agua.cgiUrl + "/agua?";

	// CREATE JSON QUERY
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "saveParameter";
	query.data = inputs;

	// SEND TO SERVER
	dojo.xhrPut(
		{
			url: url,
			contentType: "text",
			putData: dojo.toJson(query),
			timeout: 15000,
			load: function(response, ioArgs) {
				return response;
			},
			error: function(response, ioArgs) {
				return response;
			}
		}
	);

	this.savingParameter = false;

}, // Parameters.saveParameter


getFormInputs : function (widget) {
// GET INPUTS FROM THE FORM


	var inputs = new Object;

	for ( var name in this.formInputs )
	{
		var value;

		// WIDGET COMBO BOX
		if ( widget[name].value == null )
		{
			value = widget[name].getValue();
		}

		// TEXT INPUT OR HTML COMBO BOX
		else
		{
			value = widget[name].value;

			if ( widget.formInputs[name] == "word" )
			{
				value = widget.cleanWord(value);
			}
			else if ( widget.formInputs[name] == "phrase" )
			{
				value = widget.cleanEnds(value);
			}

			widget[name].value = value;
		}
		inputs[name] = value;

	}

	return inputs;
},


getEditedInputs : function (widget) {
// GET INPUTS FROM THE EDITED ITEM


	var inputs = new Object;
	for ( var name in this.formInputs )
	{

		if ( widget[name].value )
		{
			inputs[name] = String(widget[name].value.toString());

			// DON'T CLEAN SELECT BOX INPUT
			//inputs[name] = this.cleanWord(inputs[name]);
		}
		else if ( widget[name].getValue )
		{
			inputs[name] = widget[name].getValue();
		}
		else if ( widget[name].innerHTML )
		{

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

			//// REPLACE ' WITH " IN paramFunction
			//if ( name == "paramFunction" )
			//{
			//	inputs[name] = inputs[name].replace(/'/g, "\'");
			//	inputs[name] = inputs[name].replace(/[\\]+'/g, "\'");
			//}

			// INSERT THE CLEANED VALUE BACK INTO THE WIDGET
			widget[name].innerHTML = inputs[name];
		}

	}


	// SET NON-INPUT FLAG SO THESE INPUTS ARE IGNORED:
	// 	argument AND discretion
	var inputFlag = false;
	var paramType = this.paramtype.value;
	if ( paramType == 'input' )	inputFlag = true;


	// CHECK INPUTS ARE VALID AND REQUIRED INPUTS ARE NOT EMPTY
	for ( var key in this.formInputs )
	{
		// IGNORE THE argument AND discretion INPUTS IF IT'S NOT AN INPUT PARAMETER
		if ( (key == "argument" || key == "discretion")
			&& inputFlag == false )
		{
			dojo.removeClass(widget[key], 'invalid');
			continue;
		}

		var value = inputs[key];
		if ( this.isValidInput(key, value) == false
				&& this.requiredInputs[key] != null )
		{
			return null;
		}

		if ( this.requiredInputs[key] != null
				&& (value == null || value == '') )
		{

			return null;
		}
	}

	return inputs;
},


validateFormInputs : function (sourceObject, widget) {
// RETURN true IF ALL INPUTS ARE VALID, false OTHERWISE.
// ALSO ADD 'invalid' CLASS TO INVALID INPUT NODES


	// SET NON-INPUT FLAG SO THESE INPUTS ARE IGNORED:
	// 	argument AND discretion
	var inputFlag = false;

	var paramType = this.paramtype.value;
	if ( paramType == 'input' || paramType == 'Paramtype' )
		inputFlag = true;

	var validateFormInputs = true;	
	for ( var key in this.invalidInputs )
	{
		// IGNORE THE argument AND discretion INPUTS IF IT'S NOT AN INPUT PARAMETER
		if ( (key == "argument" || key == "discretion")
			&& inputFlag == false )
		{
			dojo.removeClass(widget[key], 'invalid');
			continue;
		}

		var value = sourceObject[key];

		// CHECK INPUTS ARE VALID AND REQUIRED INPUTS ARE NOT EMPTY
		if ( (widget.isValidInput(key, value) == false
				&& widget.requiredInputs[key] != null)
			|| (widget.requiredInputs[key] != null
				&& (value == null || value == '') ) )
		{
			dojo.addClass(widget[key], 'invalid');
			validateFormInputs = false;
		}
		else{
			dojo.removeClass(widget[key], 'invalid');
		}
	}

	return validateFormInputs;
},


isValidInput : function (name, value) {
	if ( this.invalidInputs[name] != null
		&& this.invalidInputs[name] == value )
	{
		return false;
	}

	return true;
},

formInputs : {
// FORM INPUTS AND TYPES (word|phrase)
	name: "word",
	argument: "word",
	type: "word",
	category: "word",
	value: "word",
	discretion: "word",
	description: "phrase",
	paramtype: "paramtype",
	format: "word",
	args: "word",
	params: "phrase",
	paramFunction: "phrase"
},


defaultInputs : {
// DEFAULT INPUTS
	name : "Name",
	argument : "Argument", 
	//type : "Type", 
	category: "Category",
	value: "Value",
	//discretion: "Discretion",
	description: "Description",
	format: "Format",
	//paramtype: "Paramtype",
	args: "Args",
	params: "Params",
	paramFunction: "ParamFunction"
},


requiredInputs : {
// REQUIRED INPUTS CANNOT BE ''
	name : 1,
	type : 1, 
	argument: 1, 
	type: 1,
	category: 1,
	discretion: 1,
	paramtype: 1
},


invalidInputs : {
// THESE INPUTS ARE INVALID
	name : "Name",
	argument : "Argument", 
	type : "Type", 
	category: "Category",
	value: "Value",
	discretion: "Discretion",
	description: "Description",
	paramtype: "Paramtype",
	format: "Format",
	args: "Args",
	params: "Params",
	paramFunction: "ParamFunction"
}


}); // plugins.admin.Parameters


/*

plug// SET LISTENERS TO ACTIVATED SAVE BUTTON AND TO CLEAR DEFAULT TEXT
// WHEN INPUTS ARE CLICKED ON
setForm : function ()
{

	// SET ADD PARAMETER ONCLICK
	dojo.connect(this.addParameterButton, "onclick", dojo.hitch(this, "saveParameter", null));

	// SET ONCLICK TO CANCEL DEFAULT TEXT
	for ( var name in this.formInputs )
	{
		dojo.connect(this[name], "onclick", dojo.hitch(this, "clearValue", this[name], this.firstLetterUpperCase(name)));
		dojo.connect(this[name], "onfocus", dojo.hitch(this, "clearValue", this[name], this.firstLetterUpperCase(name)));
	}

	//// SET NEW PARAMETER LISTENER
	//var parametersObject = this;
	//dojo.connect(this.newName, "onkeypress", function(evt){
	//	var key = evt.charOrCode;
	//	if ( key == 13 )
	//	{				
	//		// SHIFT FOCUS TO NEXT INPUT
	//		parametersObject.newDescription.focus();
	//	}
	//	
	//	// REMOVE ANY LINE RETURNS
	//	this.value = this.value.replace(/\n/, '');
	//});
	//
	//// SET NEW PARAMETER LISTENER
	//var parametersObject = this;
	//dojo.connect(this.newDescription, "onkeypress", function(evt){
	//	var key = evt.charOrCode;
	//	if ( key == 13 )
	//	{
	//		// SHIFT FOCUS TO NEXT INPUT
	//		parametersObject.newFormat.focus();
	//	}
	//
	//	// REMOVE ANY LINE RETURNS
	//	this.value = this.value.replace(/\n/, '');
	//});
	//
	//// SET NEW PARAMETER LISTENER
	//var parametersObject = this;
	//dojo.connect(this.newFormat, "onkeypress", function(evt){
	//	var key = evt.charOrCode;
	//	if ( key == 13 )
	//	{
	//		// SHIFT FOCUS TO SAVE BUTTON
	//		parametersObject.addParameterButton.focus();
	//	}
	//
	//	// REMOVE ANY LINE RETURNS
	//	this.value = this.value.replace(/\n/, '');
	//});
},



*/

