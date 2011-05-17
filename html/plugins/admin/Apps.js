dojo.provide("plugins.admin.Apps");

// ALLOW THE USER TO ADD, REMOVE AND MODIFY APPS

//dojo.require("dijit.dijit"); // optimize: load dijit layer
dojo.require("dijit.layout.ContentPane");
dojo.require("dijit.form.CheckBox");
dojo.require("dijit.form.Button");
dojo.require("dijit.form.TextBox");
dojo.require("dijit.form.Textarea");
dojo.require("dojo.parser");
dojo.require("dojo.dnd.Source");
dojo.require("plugins.core.Common");

// HAS A
dojo.require("plugins.admin.AppRow");

dojo.declare(
    "plugins.admin.Apps",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "admin/templates/apps.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

//addingApp STATE
addingApp : false,

// OR USE @import IN HTML TEMPLATE
cssFiles : [ "plugins/admin/css/apps.css" ],

// PARENT WIDGET
parentWidget : null,

	/////}}

constructor : function(args) {
	// GET INFO FROM ARGS
	this.parentWidget = args.parentWidget;
	this.apps = args.parentWidget.apps;

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
	this.tabContainer.addChild(this.appsTab);
	this.tabContainer.selectChild(this.appsTab);

	// SET APPS COMBO - WILL CASCADE TO setDragSource
	this.setAppsCombo();

	// SET NEW APP FORM
	this.setForm();

	// SET TRASH DROP TARGET
	this.setTrash();

	// SUBSCRIBE TO UPDATES
	Agua.updater.subscribe(this, "updateApps");
},


updateApps : function (args) {
// RELOAD GROUP COMBO AND DRAG SOURCE AFTER CHANGES
// TO SOURCES OR GROUPS DATA IN OTHER TABS

	// SET APPS COMBO
	this.setAppsCombo();

	// SET DRAG SOURCE
	this.setDragSource();
},

reloadAppData : function () {
// RELOAD RELEVANT DISPLAYS IN GROUP-RELATED TABS


	Agua.updater.update("updateApps");

	//var tabPaneNames = ["plugins.admin.Parameters"];
	//for ( var i in tabPaneNames )
	//{
	//	if ( this.parentWidget.paneWidgets[tabPaneNames[i]] != null )
	//	{
	//		this.parentWidget.paneWidgets[tabPaneNames[i]].reload();
	//	}
	//}
},

setAppsCombo : function (type) {
	// SET APPS COMBO BOX WITH APPLICATION NAME VALUES

	// GET APPS NAMES		
	var apps = dojo.clone(Agua.apps).custom;

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
	typesArray.splice(0,0, 'Order by Type (A-Z)');
	typesArray.splice(0,0, 'Order by Name (A-Z)');


	// SET STORE
	var data = {identifier: "name", items: []};
	for ( var i = 0; i < typesArray.length; i++ )
	{
		data.items[i] = { name: typesArray[i]	};
	}
	var store = new dojo.data.ItemFileWriteStore(	{	data: data	}	);

	this.appsCombo.popupClass = "apps appsCombo dijitReset dijitMenu";
	this.appsCombo.wrapperClass = "apps dijitPopup";
	this.appsCombo.itemHeight = 30;

	// SET COMBO
	this.appsCombo.store = store;
	this.appsCombo.startup();

	// SET COMBO VALUE
	var firstValue = typesArray[0];
	this.appsCombo.set('value', firstValue);

	//// CONNECT ONCLICK WITH dojo.connect TO BUILD TABLE
	//var appsObject = this;
	//dojo.hitch(this.appsCombo, "onChange", dojo.hitch(this, function() {
	//	this.setDragSource
	//}));

	this.setDragSource();
},

setForm : function () {
/* SET LISTENERS TO ACTIVATE SAVE BUTTON
 AND CLEAR DEFAULT TEXT WHEN INPUTS
 ARE CLICKED ON
*/


	// DEBUG
	//dojo.connect(this.newLocalOnly, "onChange", dojo.hitch(this, "submitChange", null));


	// SET ADD APP ONCLICK
	//dojo.connect(this.addAppButton, "onclick", dojo.hitch(this, "newApp", null));

	// SET ONCLICK TO CANCEL DEFAULT TEXT
	var clickableArray = [ "Name", "Type", "Location", "Executor", "Version", "Description", "Notes" ];
	for ( var i in clickableArray )
	{
		var nodeName = "new" + clickableArray[i];
		dojo.connect(this[nodeName], "onclick", dojo.hitch(this, "clearValue", this[nodeName], clickableArray[i]));
		dojo.connect(this[nodeName], "onfocus", dojo.hitch(this, "clearValue", this[nodeName], clickableArray[i]));
	}
},

setDragSource : function () {

	// DELETE EXISTING TABLE CONTENT
	while ( this.dragSource.firstChild )
	{
		this.dragSource.removeChild(this.dragSource.firstChild);
	}

	var dataArray = new Array;
	var sourceArray = dojo.clone(Agua.apps).custom;

	// FILTER SOURCE ARRAY BY type
	var type = this.appsCombo.get('value');
	if ( type == "Order by Name (A-Z)" )
	{
		sourceArray = this.sortHasharray(sourceArray, 'name');
	}
	else if ( type == "Order by Type (A-Z)" )
	{
		sourceArray = this.sortHasharray(sourceArray, 'type');
	}
	else
	{
		for ( var i = 0; i < sourceArray.length; i++ )
		{
			if ( sourceArray[i].type != type )
			{
				sourceArray.splice(i, 1);
				i--;
			}
		}

	}

	// CHECK sourceArray IS NOT NULL OR EMPTY
	if ( sourceArray == null || sourceArray.length == 0 )
	{
		return;
	}

	// GENERATE dataArray TO INSERT INTO DND APP TABLE
	for ( var j = 0; j < sourceArray.length; j++ )
	{
		var data = sourceArray[j];				
		dataArray.push( { data: data, type: ["draggableItem"] } );
	}

	// GENERATE DND APP
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


		// SET NODE name AND description
		node.name = dataArray[k].data.name;
		node.type = dataArray[k].data.type;
		node.executor = dataArray[k].data.executor;
		node.version = dataArray[k].data.version;
		node.location = dataArray[k].data.location;
		node.localonly = dataArray[k].data.localonly;
		node.description = dataArray[k].data.description;
		node.notes = dataArray[k].data.notes;
		if ( node.description == null || node.description == '' )
		{
			node.description = "Description";
		}
		if ( node.notes == null || node.notes == '' )
		{
			node.notes = "Notes";
		}

		// SET SOURCE HASH USED TO INSTANTIATE plugins.admin.AppRow
		var source = {
			name : node.name,
			type : node.type,
			executor : node.executor,
			version : node.version,
			location : node.location,
			localonly : node.localonly,
			description : node.description,
			notes : node.notes
		};
		source.parentWidget = this;

		// INSTANTIATE plugins.admin.AppRow AND APPEND TO NODE
		var sourceRow = new plugins.admin.AppRow(source);
		node.innerHTML = '';
		node.appendChild(sourceRow.domNode);
	}

	var sourceObject = this;
	dragSource.creator = function (item, hint)
	{

		var node = dojo.doc.createElement("div");
		node.name = item.name;
		node.description = item.description;
		node.notes = item.notes;
		node.id = dojo.dnd.getUniqueId();
		node.className = "dojoDndItem";


		// SET FANCY FORMAT IN NODE INNERHTML
		node.innerHTML = "<table> <tr><td><strong style='color: darkred'>" + item.name + "</strong></td></tr><tr><td> " + item.description + "</td></tr></table>";

		return {node: node, data: item, type: ["text"]};
	};

},
editAppRow : function (appRowWidget, node) {
/* EDIT APPLICATION ROW - SAVE IF 'RETURN' PRESSED,
 EXIT WITHOUT CHANGES IF 'ESCAPE' PRESSED */


	//////console.dir(node);


	// RETURN IF ALREADY EDITING APP ROW (I.E., MULTIPLE CLICKS)
	if ( this.editingAppRow == true ) return;
	this.editingAppRow = true;

	// REPLACE THE TD INNERHTML WITH A TEXTAREA
	var text = node.innerHTML;
	if ( text.match(/App notes/) )
	{
		text = '';
	}

	node.innerHTML = '';
	if ( text == null || ! text ) text = '';

	// RETURN IF THIS IS A DOUBLE-CLICK
	if ( text.match(/^<i/) ||
		text.match(/^<br/) ||
		text.match(/^<fieldset/) ||
		text.match(/^<textarea/) )
	{
		this.editingAppRow = false;
		return;
	}

	// CREATE INPUT TEXT AREA
	var textarea = document.createElement('textarea');
	dojo.addClass(textarea, 'editAppRow');
	node.appendChild(textarea);
	textarea.value = text;
	textarea.focus();

	// SET NEW APP LISTENER
	var appsObject = this;
	dojo.connect(textarea, "onkeypress", function(evt){

		// summary: handles keyboard events
		var key = evt.charOrCode;
		evt.stopPropagation();

		if ( key == 13 )
		{
			var value = textarea.value;

			// REMOVE TEXTAREA
			node.removeChild(textarea);
			node.innerHTML = value;

			// NOT EDITING ROW ANY MORE
			appsObject.editingAppRow = false;

			// IF INPUT IS INVALID, ADD 'invalid' CSS CLASS AND QUIT
			var key = node.getAttribute('class').match(/^(\S+)/)[1];
			if ( (appsObject.isValidInput(key, value) == false
					&& appsObject.requiredInputs[key] != null)
				|| (appsObject.requiredInputs[key] != null
					&& value == null || value == '') )
			{
				dojo.addClass(node, 'invalid');

				// RESTORE ORIGINAL VALUE OF INPUT
				node.innerHTML = text;
			}

			// OTHERWISE, SAVE THE APP AND RELOAD THE DRAG SOURCE
			else
			{
				dojo.removeClass(node, 'invalid');

				// GET INPUTS
				var inputs = appsObject.getEditedInputs(appRowWidget);
				if ( inputs == null ) return;

				// SAVE APPLICATION
				appsObject.saveApp(inputs);
			}
		}
		else if (key == dojo.keys.ESCAPE)
		{
			appsObject.editingAppRow = false;

			// REMOVE TEXTAREA
			node.removeChild(textarea);

			// RESTORE ORIGINAL VALUE
			node.innerHTML = text;
		}
	});

	// QUIT EDIT IF FOCUS IS LOST
	dojo.connect(appRowWidget, "onBlur", function(evt){
		appsObject.editingParameterRow = false;

		// REMOVE TEXTAREA
		node.removeChild(textarea);

		// RESTORE ORIGINAL VALUE
		node.innerHTML = text;
	});

	// CLEAN UP WHITESPACE
	appObject.name = appObject.name.replace(/\s+$/,'');
	appObject.name = appObject.name.replace(/^\s+/,'');

	// REMOVING APP FROM Agua.apps
	Agua.removeApp(appObject, "custom");

	// RESET THE APPS TABLE
	this.setDragSource();

	var url = Agua.cgiUrl + "/agua?";

	// CREATE JSON QUERY
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "deleteApp";
	query.data = appObject;

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

	// RELOAD RELEVANT DISPLAYS IN GROUP-RELATED TABS
	Agua.updater.update("updateApps");

}, // Apps.deleteApp


deleteApp : function (appObject) {
// DELETE APPLICATION FROM Agua.apps OBJECT AND IN REMOTE DATABASE

	// CLEAN UP WHITESPACE
	appObject.name = appObject.name.replace(/\s+$/,'');
	appObject.name = appObject.name.replace(/^\s+/,'');

	// REMOVING APP FROM Agua.apps
	Agua.removeApp(appObject, "custom");

	// RESET THE APPS TABLE
	this.setDragSource();

	var url = Agua.cgiUrl + "/agua?";

	// CREATE JSON QUERY
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "deleteApp";
	query.data = appObject;

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

	// RELOAD RELEVANT DISPLAYS IN GROUP-RELATED TABS
	Agua.updater.update("updateApps");

}, // Apps.deleteApp


addApp : function (sourceObject) {
//	ADD APPLICATION TO Agua.apps AND TO REMOTE DATABASE
//return;

	if ( Agua.isApp(appObject.name) )
	{
		dojo.removeClass(this.newName, 'invalid');
		dojo.addClass(this.newName, 'invalid');
		this.newName.focus();
		return;
	}

	// ADD APP TO Agua.apps ARRAY
	Agua.addApp(appObject, "custom");

	// RESET APP TABLE
	this.setDragSource();

	// SAVE APP TO REMOTE DATABASE
	this.newApp(appObject);

	// RELOAD RELEVANT DISPLAYS IN GROUP-RELATED TABS
	this.reloadAppData();

},	// Apps.addApp


saveApp : function (inputs) {
//	SAVE AN APPLICATION TO Agua.apps AND TO REMOTE DATABASE

	if ( this.savingApp == true )	return;
	this.savingApp = true;

	if ( inputs == null )
	{
		inputs = this.getFormInputs();

		// RETURN IF INPUTS ARE NULL OR INVALID
		if ( inputs == null || this.validateFormInputs(inputs) == false )
		{
			this.savingApp = false;
			return;
		}
	}


	// OTHERWISE, REMOVE ORIGINAL APPLICATION OBJECT FROM Agua.apps 
	// THEN ADD NEW APPLICATION OBJECT TO Agua.apps ARRAY
	Agua.removeApp({ name: inputs.name } , "custom");
	Agua.addApp(inputs, "custom");

	// REDO APP TABLE
	//this.setDragSource();

	// SAVE NEW APP TO REMOTE DATABASE
	var url = Agua.cgiUrl + "/agua?";

	// CREATE JSON QUERY
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "saveApp";
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

	this.savingApp = false;

	// RELOAD RELEVANT DISPLAYS IN GROUP-RELATED TABS
	this.reloadAppData();

}, // Apps.saveApp


newApp : function (event) {
//	SAVE AN APPLICATION TO Agua.apps AND TO REMOTE DATABASE

	if ( this.savingApp == true )	return;
	this.savingApp = true;

	var inputs = this.getFormInputs();

	// RETURN IF INPUTS ARE NULL OR INVALID
	if ( inputs == null || this.validateFormInputs(inputs) == false )
	{
		this.savingApp = false;
		return;
	}


	// OTHERWISE, REMOVE ORIGINAL APPLICATION OBJECT FROM Agua.apps 
	// THEN ADD NEW APPLICATION OBJECT TO Agua.apps ARRAY
	Agua.removeApp({ name: inputs.name }, "custom" );
	Agua.addApp(inputs, "custom");

	// REDO APP TABLE
	this.setDragSource();

	// CREATE APP OBJECT
	var appObject = new Object;
	appObject.name = inputs.name;
	appObject.type = inputs.type;
	appObject.executor = inputs.executor;
	appObject.version = inputs.version;
	appObject.location = inputs.location;
	appObject.localonly = inputs.localonly;
	appObject.description = inputs.description;
	appObject.notes = inputs.notes;


	// SAVE NEW APP TO REMOTE DATABASE
	var url = Agua.cgiUrl + "/agua?";

	// CREATE JSON QUERY
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "saveApp";
	query.data = appObject;

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

	this.savingApp = false;

	// RELOAD RELEVANT DISPLAYS IN GROUP-RELATED TABS
	this.reloadAppData();

	// RESET THE APPS TABLE
	this.setDragSource();

}, // Apps.newApp



setTrash : function () {
//	DELETE NODE IF DROPPED INTO TRASH. ACTUAL REMOVAL FROM THE
//	DATA IS ACCOMPLISHED IN THE onDndDrop LISTENER OF THE APP

	var trash = new dojo.dnd.Source(
		this.trashContainer,
		{
			accept : [ "draggableItem" ]
		}
	);

	// REMOVE DUPLICATE NODES
	var appsObject = this;
	dojo.connect(trash, "onDndDrop", function(source, nodes, copy, target){


		// NODE DROPPED ON SELF --> DELETE THE NODE
		if ( target == this )
		{
			for ( var i = 0; i < nodes.length; i++ )
			{
				var node = nodes[i];

				// HACK TO AVOID THIS ERROR: node.parentNode is null
				try {


					appsObject.deleteApp({
						name: node.name,
						type: node.type,
						location: node.location
					});

				}
				catch (e)
				{
				}
			}

			// DELETE EXISTING TABLE CONTENT
			////console.dir(appsObject.trashContainer);
			while ( appsObject.trashContainer.childNodes.length > 2 )
			{

				appsObject.trashContainer.removeChild(appsObject.trashContainer.childNodes[2]);
			}
		}
	});

}, // setTrash


clearValue : function (textarea, value) {
	if ( textarea.value == value )
	{
		textarea.value = '';
	}
},


getFormInputs : function () {
// GET INPUTS FROM THE FORM
	////console.clear();

	var inputs = new Object;
	var inputsArray = [ "name", "type", "executor", "version", "location", "localonly", "description", "notes" ];
	for ( var i in inputsArray )
	{
		var nodeName = "new" + this.firstLetterUpperCase(inputsArray[i]);
		var value;

		if ( inputsArray[i] == "localonly" )
		{
			if ( this[nodeName].getValue() == "on" )
			{
				value = 1;
			}
			else
			{
				value = 0;
			}

		}
		else
		{
			value = this[nodeName].value;
		}

		inputs[inputsArray[i]] = value;
	}

	// REMOVE WHITESPACE
	inputs.name = this.cleanWord(inputs.name);
	inputs.type = this.cleanWord(inputs.type);
	inputs.version = this.cleanWord(inputs.version);
	inputs.location = this.cleanWord(inputs.location);

	inputs.executor = this.cleanEnds(inputs.executor);
	inputs.description = this.cleanEnds(inputs.description);
	inputs.notes = this.cleanEnds(inputs.notes);

	return inputs;
},



getEditedInputs : function (widget) {
// GET INPUTS FROM THE EDITED ITEM

	var inputs = new Object;
	var inputsArray = [ "name", "type", "executor", "version", "location", "localonly", "description", "notes" ];
	for ( var i in inputsArray )
	{
		if ( inputsArray[i] == "localonly" )
		{
			inputs[inputsArray[i]] = 0;
			if ( widget[inputsArray[i]].getValue() == "on" )
			{
				inputs[inputsArray[i]] = 1;
			}
		}
		else
		{
			inputs[inputsArray[i]] = widget[inputsArray[i]].innerHTML;
		}
	}

	// REMOVE WHITESPACE
	inputs.name = this.cleanWord(inputs.name);
	inputs.type = this.cleanWord(inputs.type);
	inputs.version = this.cleanWord(inputs.version);
	inputs.location = this.cleanWord(inputs.location);

	inputs.executor = this.cleanEnds(inputs.executor);
	inputs.description = this.cleanEnds(inputs.description);
	inputs.notes = this.cleanEnds(inputs.notes);


	// CHECK INPUTS ARE VALID AND REQUIRED INPUTS ARE NOT EMPTY
	for ( var i in inputsArray )
	{
		var key = inputsArray[i];
		var value = inputs[inputsArray[i]];
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



validateFormInputs : function (sourceObject, node) {
// RETURN true IF ALL INPUTS ARE VALID, false OTHERWISE.
// ALSO ADD 'invalid' CLASS TO INVALID INPUT NODES


	var validateFormInputs = true;
	var inputsArray = [ "name", "type", "executor", "version", "location", "description", "notes" ];
	var invalid = [ "Name", "Type", "Executor", "Location"];

	for ( var key in sourceObject )
	{
		var nodeName = "new" + this.firstLetterUpperCase(key);
		var value = sourceObject[key];

		// CHECK INPUTS ARE VALID AND REQUIRED INPUTS ARE NOT EMPTY
		if ( (this.isValidInput(key, value) == false
				&& this.requiredInputs[key] != null)
			|| (this.requiredInputs[key] != null
				&& (value == null || value == '') ) )
		{
			dojo.addClass(this[nodeName], 'invalid');
			validateFormInputs = false;
		}
		else{
			dojo.removeClass(this[nodeName], 'invalid');
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

requiredInputs : {
// REQUIRED INPUTS CANNOT BE ''
	name : 1,
	type : 1, 
	location: 1
},

invalidInputs : {
// THESE INPUTS ARE INVALID
	name : "Name",
	type : "Type", 
	executor: "Executor", 
	version: "Version",
	location: "Location",
	description: "Description",
	notes: "Notes"
},

defaultInputs : {
// THESE INPUTS ARE default
	name : "Name",
	type : "Type", 
	executor: "Executor", 
	version: "Version",
	location: "Location",
	description: "Description",
	notes: "Notes"
}



}); // plugins.admin.Apps

