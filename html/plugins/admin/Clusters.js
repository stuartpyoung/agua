dojo.provide("plugins.admin.Clusters");

// ALLOW THE USER TO ADD, REMOVE AND MODIFY StarCluster CLUSTER GROUPS

//dojo.require("dijit.dijit"); // optimize: load dijit layer
dojo.require("dijit.form.Button");
dojo.require("dijit.form.TextBox");
dojo.require("dijit.form.Textarea");
dojo.require("dojo.parser");
dojo.require("dojo.dnd.Source");
dojo.require("plugins.core.Common");

// FORM VALIDATION
dojo.require("plugins.form.ValidationTextarea");

// HAS A
dojo.require("plugins.admin.ClusterRow");

dojo.declare(
    "plugins.admin.Clusters",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	/////}	

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "admin/templates/clusters.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

//addingUser STATE
addingUser : false,

// OR USE @import IN HTML TEMPLATE
cssFiles : [ dojo.moduleUrl("plugins") + "/admin/css/clusters.css"],

// PARENT WIDGET
parentWidget : null,

// FORM INPUTS AND TYPES (word|phrase)
formInputs : {
	clusterName	:	"word",
	minNodes	:	"number",
	maxNodes	:	"number",
	typeCombo	:	"combo",
	amiId		:	"word",
	description	:	"phrase",
	notes		:	"phrase"
},

defaultInputs : {
// DEFAULT INPUTS
	clusterName	:	"clusterName",
	minNodes	:	"minNodes",
	maxNodes	:	"maxNodes",
	instanceType:	"instanceType", 
	amiId		:	"AMI ID",
	description	:	"Description",
	notes		:	"Notes"
},

requiredInputs : {
// REQUIRED INPUTS CANNOT BE ''
// combo INPUTS ARE AUTOMATICALLY NOT ''
	clusterName : 1,
	minNodes 	: 1, 
	maxNodes	: 1, 
	instanceType: 1,
	amiId		: 1
},

invalidInputs : {
// THESE INPUTS ARE INVALID
	clusterName : 	"myCluster",
	minNodes	:	"minNodes",
	maxNodes	:	"0",
	amiId		:	"AMI ID",
	description	:	"Description",
	notes		:	"Notes"
},

// STARTUP METHODS
constructor : function(args) {
	// GET INFO FROM ARGS
	this.parentWidget = args.parentWidget;
	this.clusters = args.parentWidget.clusters;

	// LOAD SORIA AND FILEPICKER CSS
	this.loadCSS();		
},

postCreate : function() {

	this.startup();
},

startup : function () {

	// COMPLETE CONSTRUCTION OF OBJECT
	this.inherited(arguments);	 

	// ADD ADMIN TAB TO TAB CONTAINER		
	this.tabContainer.addChild(this.clustersTab);
	this.tabContainer.selectChild(this.clustersTab);

	// SET DRAG SOURCE - LIST OF CLUSTERS

	// SET NEW PARAMETER FORM
	this.setForm();

	this.setDragSource();

	// SET TRASH DROP TARGET
	this.setTrash();

	// SUBSCRIBE TO UPDATES
	Agua.updater.subscribe(this, "updateClusters");
},

updateClusters : function (args) {
// RELOAD THE GROUP COMBO AND DRAG SOURCE AFTER CHANGES
// TO SOURCES OR GROUPS DATA IN OTHER TABS


	// SET DRAG SOURCE
	this.setDragSource();
},


setForm : function () {
// SET LISTENERS TO ACTIVATED SAVE BUTTON AND TO CLEAR DEFAULT TEXT
// WHEN INPUTS ARE CLICKED ON

	// SET ADD PARAMETER ONCLICK
	dojo.connect(this.addParameterButton, "onclick", dojo.hitch(this, "saveParameter", null));	
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

				//console.dir(node);

				// HACK TO AVOID THIS ERROR: node.parentNode is null
				try {
					node.parentNode.removeChild(node);


					var sourceObject = new Object;
					sourceObject.cluster = node.clusterName;
					sourceObject.instancetype = node.instanceType;
					sourceObject.amiid = node.amiId;
					sourceObject.minnodes = node.minNodes;
					sourceObject.maxnodes = node.maxNodes;
					sourceObject.description = node.description;
					sourceObject.notes = node.notes;

					thisObject.removeCluster(sourceObject);
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

setDragSource : function () {
// SET THE DRAG SOURCE WITH PARAMETER OBJECTS

	// DELETE EXISTING CONTENT
	while ( this.dragSource.firstChild )
	{
		this.dragSource.removeChild(this.dragSource.firstChild);
	}

	// INITIALISE USER INFO
	var clusters = Agua.getClusters();
	clusters = this.sortHasharray(clusters, 'cluster');

	// CHECK IF clusters IS NULL
	if ( clusters == null )
	{
		return;
	}

	// GENERATE dataArray TO INSERT INTO DND PARAMETER TABLE
	var dataArray = new Array;
	for ( var i in clusters )
	{
		var data = clusters[i];				
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
		// SET SOURCE HASH USED TO INSTANTIATE plugins.admin.ClusterRow
		var nameArray = [ "cluster", "minnodes", "maxnodes", "instancetype", "amiid", "description", "notes" ];
		var valueArray = [ "clusterName", "minNodes", "maxNodes", "instanceType", "amiId", "description", "notes" ];
		var dataObject = new Object;
		for ( var i in valueArray )
		{
			node[valueArray[i]] = dataObject[valueArray[i]] = dataArray[k].data[nameArray[i]] ? dataArray[k].data[nameArray[i]] : this.firstLetterUpperCase(valueArray[i]);
		}
		dataObject.parentWidget = this;

		// INSTANTIATE plugins.admin.ClusterRow
		var clusterRow = new plugins.admin.ClusterRow(dataObject);
		node.innerHTML = '';
		node.appendChild(clusterRow.domNode);
		node.widget = clusterRow;
	}


	var sourceObject = this;
	dragSource.creator = function (item, hint)
	{

		var node = dojo.doc.createElement("div");
		node.name = item.cluster;
		node.id = dojo.dnd.getUniqueId();
		node.className = "dojoDndItem";


		// SET FANCY FORMAT IN NODE INNERHTML
		node.innerHTML = "<table> <tr><td><strong style='color: darkred'>" + item.cluster + "</strong></td></tr><tr><td> " + item.description + "</td></tr></table>";

		return {node: node, data: item, type: ["text"]};
	};
},

// ADD/REMOVE CLUSTER METHODS
saveCluster : function (inputs, reload) {
//	SAVE A PARAMETER TO Agua.parameters AND TO REMOTE DATABASE


	if ( this.savingCluster == true )	return;
	this.savingCluster = true;

	if ( inputs == null )
	{
		inputs = this.getFormInputs(this);

		// RETURN IF INPUTS ARE NULL OR INVALID
		if ( inputs == null || this.validateFormInputs(inputs, this) == false )
		{
			this.savingCluster = false;
			return;
		}
	}

	var clusterObject = new Object;
	clusterObject.username = Agua.cookie('username');	
	clusterObject.sessionId = Agua.cookie('sessionId');	
	clusterObject.cluster = inputs.clusterName;
	for ( var input in inputs )
	{
		clusterObject[input.toLowerCase()] = inputs[input];
	}
	clusterObject.instancetype = inputs["typeCombo"];

	Agua.removeCluster(clusterObject); // not necessary as addCluster does delete first
	Agua.addCluster(clusterObject);

	this.savingCluster = false;

	// REDO PARAMETER TABLE
	if ( reload != null )	this.setDragSource();

	// SAVE ON REMOTE DATABASE
	var url = Agua.cgiUrl + "/agua?";
	clusterObject.mode = "addCluster";

	// SEND TO SERVER
	dojo.xhrPut(
		{
			url: url,
			contentType: "text",
			putData: dojo.toJson(clusterObject),
			timeout: 15000,
			load: function(response, ioArgs) {
				return response;
			},
			error: function(response, ioArgs) {
				return response;
			}
		}
	);


}, // Clusters.saveCluster


addCluster : function (event) {

	if ( this.savingCluster == true )
	{
		return;
	}
	this.savingCluster = true;

	var inputs = this.getFormInputs(this);

	// RETURN IF INPUTS ARE NULL OR INVALID
	if ( inputs == null || this.validateFormInputs(inputs, this) == false )
	{
		this.savingCluster = false;
		return;
	}

	// GENERATE CLUSTER OBJECT
	var clusterObject = new Object;
	clusterObject.username = Agua.cookie('username');	
	clusterObject.sessionId = Agua.cookie('sessionId');	
	clusterObject.mode = "addCluster";	
	clusterObject.cluster = inputs.clusterName;
	for ( var input in inputs )
	{
		clusterObject[input.toLowerCase()] = inputs[input];
	}
	clusterObject.instancetype = inputs["typeCombo"];

	// ADD TO Agua LOCALLY
	Agua.addCluster(clusterObject);

	// SEND TO SERVER
	var url = Agua.cgiUrl + "/agua?";
	dojo.xhrPut(
		{
			url: url,
			contentType: "text",
			putData: dojo.toJson(clusterObject),
			timeout: 15000,
			load: function(response, ioArgs) {
				return response;
			},
			error: function(response, ioArgs) {
				return response;
			}
		}
	);

	this.savingCluster = false;

	// RESET PARAMETER TABLE
	this.setDragSource();

}, // Clusters.addCluster

removeCluster : function (clusterObject) {
// DELETE PARAMETER FROM Agua.parameters OBJECT AND IN REMOTE DATABASE


	// REMOVING PARAMETER FROM Agua.parameters
	var success = Agua.removeCluster(clusterObject)
	if ( success == false )
	{
		return;
	}

	// RESET THE PARAMETERS TABLE
	this.setDragSource();

	var url = Agua.cgiUrl + "/agua?";

	// CREATE JSON QUERY
	clusterObject.username = Agua.cookie('username');
	clusterObject.sessionId = Agua.cookie('sessionId');
	clusterObject.mode = "removeCluster";

	// SEND TO SERVER
	dojo.xhrPut(
		{
			url: url,
			contentType: "text",
			putData: dojo.toJson(clusterObject),
			timeout: 15000,
			load: function(response, ioArgs) {
				return response;
			},
			error: function(response, ioArgs) {
				return response;
			}
		}
	);

}, // Clusters.removeCluster


// CLUSTER INPUT METHODS
editClusterRow : function (parameterRowWidget, node) {
// EDIT PARAMETER ROW - SAVE IF 'RETURN' PRESSED,
// EXIT WITHOUT CHANGES IF 'ESCAPE' PRESSED


	// RETURN IF ALREADY EDITING PARAMETER ROW (I.E., MULTIPLE CLICKS)
	if ( this.editingClusterRow == true ) return;
	this.editingClusterRow = true;

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
	dojo.addClass(textarea, 'editClusterRow');
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
			appsObject.editingClusterRow = false;

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
				appsObject.saveCluster(inputs, "reload");
			}
		}
		else if (key == dojo.keys.ESCAPE)
		{
			appsObject.editingClusterRow = false;

			// REMOVE TEXTAREA
			node.removeChild(textarea);

			// RESTORE ORIGINAL VALUE
			node.innerHTML = text;
		}
	});

	// QUIT EDIT IF FOCUS IS LOST
	dojo.connect(parameterRowWidget, "onBlur", function(evt){
		appsObject.editingClusterRow = false;

		// REMOVE TEXTAREA IF STILL ATTACHED
		if ( textarea != null && textarea.parentNode != null )
		{
			textarea.parentNode.removeChild(textarea);
		}

		// RESTORE ORIGINAL VALUE
		node.innerHTML = text;
	});
},

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

			if ( widget.formInputs[name] == "number" )
			{
				value = widget.cleanNumber(value);
			}
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
// GET INPUTS FROM THE EDITED ITEM.
// EDITED ITEM = CLUSTER FORM OR CLUSTER ROW

	//console.dir(widget);

	var inputs = new Object;
	for ( var name in this.formInputs )
	{

		if ( widget[name].id.match(/^dijit_form_NumberTextBox/) )
		{
			inputs[name] = String(widget[name]);
		}

		else if ( widget[name].get && widget[name].get('value') )
		{
			inputs[name] = String(widget[name].get('value').toString());

			// DON'T CLEAN SELECT BOX INPUT
			//inputs[name] = this.cleanWord(inputs[name]);
		}
		else if ( widget[name].value )
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
		else { 
		}

		if ( inputs[name] == null )
		{
			return;	
		}

	}

	//// SET NON-INPUT FLAG SO THESE INPUTS ARE IGNORED:
	//// 	argument AND discretion
	//var inputFlag = false;
	//var paramType = this.paramtype.value;
	//if ( paramType == 'input' )	inputFlag = true;


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


	var validateFormInputs = true;	
	for ( var key in this.invalidInputs )
	{
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

// UTILITY METHODS
cleanNumber : function (number) {

	if ( number == null )	return '';

	number = number.toString().replace(/[\s\n]+/, '');

	return number.toString().replace(/[\s\n]+/, '');
},

cleanWord : function (word) {

	if ( word == null )	return '';

	return word.replace(/[\s\n]+/, '');
},

cleanEnds : function (words) {
	if ( words == null )	return '';

	words = words.replace(/^[\s\n]+/, '');
	words = words.replace(/[\s\n]+$/, '');

	return words;
},


cleanEdges : function (string) {
// REMOVE WHITESPACE FROM EDGES OF TEXT

	if ( string == null )	{ 	return null; }
	string = string.replace(/^\s+/, '');
	string = string.replace(/\s+$/, '');
	return string;
},

checkEnter : function (event) {

	if (event.keyCode == dojo.keys.ENTER)
	{
		this.saveCluster();
		dojo.stopEvent(event);
	}
},

checkEnterNodes : function (event) {

	if (event.keyCode == dojo.keys.ENTER)
	{
		document.body.focus();

		this.checkNodeNumbers();

		this.saveCluster();
		dojo.stopEvent(event);
	}
},


checkNodeNumbers : function () {
// SET MIN NODES VALUE TO SENSIBLE NUMBER 


	if (this.minNodes.value > this.maxNodes.value )
	{
		this.minNodes.set('value', this.maxNodes.value);
	}
},

clearValue : function (event) {
	var textarea = event.target;

	if ( textarea == null )	return;
	if ( textarea.value == null )	return;

	if ( textarea.value == textarea.defaultValue )
	{
		textarea.value = '';
	}
}




}); // plugins.admin.Clusters

/*
initialiseClusters : function () {

	// DISPLAY USERNAME
	var username = Agua.cookie('username');
	this.username.innerHTML = username;

	// INITIALISE USER INFO
	var clusters = Agua.getClusters(username);
	if ( ! clusters.length > 0 )	return;
	this.minnodes.value = clusters[1].minnodes;
	this.maxnodes.value = clusters[1].maxnodes;
	this.amiId.value = clusters[1].amiid;
},

*/
