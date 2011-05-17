dojo.provide("plugins.report.SNP");

// TO DO: ADDITIONAL FILTER WITH CCDS
//
//	Promoter/Upstream by  ___  bases 
//	5' UTR Exons 
//	CDS Exons 
//	3' UTR Exons 
//	Introns 
//	Downstream by ___ bases

////http://localhost/agua/0.4/dojo-1.5.0/dojox/layout/tests/test_ExpandoPane.html
////http://localhost/agua/0.4/dojo-1.5.0/dojox/layout/tests/test_ExpandoPane.html


dojo.require("dijit.dijit"); // optimize: load dijit layer

dojo.require("dijit.form.CheckBox");
////dojo.require("dijit.form.ValidationTextBox");
//dojo.require("dijit.form.RadioButton");
//dojo.require("dijit.form.ComboBox");

dojo.require("dijit.form.Slider");
dojo.require("dijit.form.FilteringSelect");
dojo.require("dijit.form.Button");

// DEPTH NUMBER SPINNER
dojo.require("dijit.form.NumberSpinner");

// NOTES EDITOR
dojo.require("dijit.Editor");
dojo.require("dijit.form.DateTextBox");
dojo.require("dijit.form.Textarea");

dojo.require("dijit.form.TextBox");
dojo.require("dijit.form.ValidationTextBox");
dojo.require("dijit.form.NumberTextBox");
dojo.require("dijit.form.CurrencyTextBox");
dojo.require("dojo.currency");
dojo.require("dijit.Dialog");

// FILE UPLOAD
//dojo.require("plugins.upload.FileUp");
dojo.require("plugins.upload.FileUpload");

// TITLE PANES FOR EACH SECTION
dojo.require("dijit.TitlePane");

// STORE FOR PROJECT AND WORKFLOW COMBOS
dojo.require("dojo.data.ItemFileReadStore");

// GRID BUT ARE NEEDED FOR SNP ??
dojo.require("dijit.Tooltip");
dojo.require("dijit.Menu");
dojo.require("dijit.ColorPalette");

// PARSER
dojo.require("dojo.parser");

// CUSTOM EDITOR
dojo.require("plugins.report.Editor");

// FILE MANAGER HAS FILE SELECTORS
dojo.require("plugins.workflow.FileManager");
dojo.require("plugins.workflow.FileSelector");

// INHERITED CLASSES
dojo.require("plugins.core.Common");

// HAS A
dojo.require("plugins.report.Grid");
dojo.require("plugins.core.ComboBox");

dojo.declare(
    "plugins.report.SNP",
    [ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	/////}}

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "report/templates/snp.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

url : '',
	id : '',
filename: '',
loading: null,

// OR USE @import IN HTML TEMPLATE
cssFiles : [ "plugins/report/css/snp.css" ],

// CURRENT PROJECT AND WORKFLOW
project : null,
workflow : null,	  

fileBrowser : null,

// FILTERS
filters: [ "chromosome","variant","depth","sense","exonic","dbsnp" ],
filterInputs: [ "chromosomeCombo","variantInput","depthSpinner","senseCombo","exonicCombo","dbsnpCombo" ],

// ELEMENTS
elementTypes :
{
	"editor" : "editor",
	"reportCombo": "combobox",
	"workflowCombo": "combobox",
	"projectCombo": "combobox",
	"chromosomeCombo": "combobox",
	"speciesCombo": "combobox",
	"senseCombo": "combobox",
	"exonicCombo": "combobox",
	"dbsnpCombo": "combobox",
	"chromosomeCheckbox": "checkbox",
	"senseCheckbox": "checkbox",
	"exonicCheckbox": "checkbox",
	"dbsnpCheckbox": "checkbox",
	"depthCheckbox": "checkbox", 
	"variantCheckbox": "checkbox",
	"fileInput": "input",
	"variantInput": "input",
	"depthSpinner": "spinner",
	"totalResult": "div",
	"chromosomeResult": "div",
	"variantResult": "div",
	"depthResult": "div",
	"senseResult": "div",
	"exonicResult": "div",
	"dbsnpResult": "div"
},


constructor : function (args) {

	this.project = args.project;
	this.workflow = args.workflow;
	this.filename = args.filename;

},

postMixInProperties: function() {
	//this.popup = new dijit.Dialog({});
},

postCreate: function() {

	this.startup();
},

// START UP GRID OBJECT, SET COMBOS, BUTTONS, ETC.
startup : function () {

	// SET UP THE ELEMENT OBJECTS AND THEIR VALUE FUNCTIONS
	this.inherited(arguments);

	// LOAD CSS STYLE FILES
	this.loadCSS();

	// ADD TO TAB CONTAINER		
	this.attachNode.appendChild(this.containerNode);

//return;

	// SET GRID
	this.setGrid();

	// SET REPORT COMBO KEY PRESS LISTENER
	this.setReportComboOnkeypress();

	// SET PROJECT COMBO
	this.setProjectCombo();

	// SET FILE NAME
	this.setFilename();

	// SET EDITOR
	this.setEditor();

	// SET SPECIES COMBO
	this.setSpeciesCombo();

	// SET VARIANT SLIDER
	this.setVariantSlider();

	// SET DEPTH SPINNER
	this.setDepthSpinner();

	// SET COMBO BOXES
	this.setCombos();

	//// SET FILE UPLOAD
	//this.setFileUpload();

	//// SET FILE BROWSE
	//this.setFileBrowse();

	// SET 'OK' BUTTON
	this.setOkButton();	

	// SET CHECKBOXES
	this.setCheckboxes();


// DEBUG
this.filterReport();


},

// INSTANTIATE GRID OBJECT AND ATTACH
setGrid : function () {
	this.grid = new plugins.report.Grid(
		{
			//attachNode: this.gridAttachNode,
			attachNode: this.gridPane.domNode,
			parentWidget: this
		}
	);
},

// LOAD FILTER RESULT DATA INTO GRID OBJECT
loadGrid : function(data) {

	this.grid.load(data);
},


// FILTER SNPS BASED ON FILTER PARAMETERS
filterReport : function() {
	var snpReport = this;

//			// SET LOADING ICON
//            this.loading = true;
//            dojo.removeClass(this.loadingIcon, "static");
//            dojo.addClass(this.loadingIcon, "loading");

	// GENERATE QUERY
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');

	if ( query.sessionId == null || query.username == null )
	{
		return;
	}
	query.project = this.projectCombo.getValue();
	query.workflow = this.workflowCombo.getValue();
	query.report = this.reportCombo.getValue();

	// ADD mode AND class
	query.mode = "filterReport";
	query["class"] = "Report::SNP";

	// GET THE CURRENT SETTINGS OF ALL THE REPORT FORM WIDGETS
	for ( name in this.elementTypes )
	{
		var value = this.getValueType(this[name], this.elementTypes[name]);

		if ( name == "fileInput" )
		{
			var project = this.projectCombo.getValue();
			var workflow = this.workflowCombo.getValue();
			value = project + "/" + workflow + "/" + value;
		}

		if ( value == "on" )	value = true;
		if ( value == "off" )	value = false;
		if ( ! value )	value = false;
		query[name] = value;
	}


	//// GET THE FILTER VALUES
	//var filters = this.filters;
	//query.filters = new Object;
	//for ( var i = 0; i < filters.length; i++ )
	//{
	//	var name = filters[i] + "Checkbox";
	//	var value;
	//	var value = this.getValueType(this[filters[i]], this.elementTypes[filters[i]]);
	//	if ( value == "on" )	{	value = true;	}
	//	if ( value == "off" )	{	value = false;	}
	//	query.filters[filters[i]] = value;
	//}


	// POST DATA TO SERVER
	var template = this;
	dojo.xhrPut(
		{
			//url: this.url,
			url : Agua.cgiUrl + "report.cgi",
			contentType: "text",
			handleAs: "json",
			sync: true,
			postData: dojo.toJson(query),
			//handleAs: "json-comment-filtered",
			timeout: 50000,
			load: function(response, ioArgs) {
				snpReport.loadReturnValues(response);
			},
			error: function(response, ioArgs) {
				return response;
			}
		}
	);
},


// LOAD RETURNED VALUES FROM SERVER INTO WIDGETS/NODES ON PAGE
loadReturnValues : function (values) {

	// SET UP FILTERS RESULTS IN ORDER
	for ( var i = 0; i < this.filters.length; i++)
	{
		var filterName = this.filters[i] + "Result";
		var value = values[filterName];
		if ( value == null )	value = "0";
		this.setValueType(this[filterName], this.elementTypes[filterName], value);
	}

	// DO totalResult
	var totalValue = values["totalResult"] ? values["totalResult"] : "0";
	this.setValueType(this.totalResult, this.elementTypes["totalResult"], totalValue);

	// DO GRID
	var data = {"identifier":"id","label":"id","items":[]};

	// CONVERT ARRAYS INTO HASHES
	var fieldNames = ['name','chromosome','ccdsstart','ccdsstop','referencenucleotide','variantnucleotide','depth','variantfrequency','chromosomestart','chromosomestop','sense','referencecodon','variantcodon','referenceaa','variantaa','strand','snp','score','dbsnpstrand'];


	// SET COMBOBOX VALUES
	var twoDarray = dojo.fromJson(values["outputResult"]);
	if ( twoDarray == null )
	{
		return;
	}
	for ( var j = 0; j < twoDarray.length; j++ )
	{
		var hash = new Object;
		hash.id = j + 1;
		var array = twoDarray[j];

		for ( k = 0; k < fieldNames.length; k++ )
		{
			hash[fieldNames[k]] = array[k];
		}
		data.items.push(hash);
	}


	// LOAD GRID
	setTimeout(function(thisObj) { thisObj.loadGrid(data); }, 100, this);

	//// SET FINISHED LOADING ICON
	//this.loading = false;
	//dojo.removeClass(this.loadingIcon, "loading");
	//dojo.addClass(this.loadingIcon, "static");
},


// SET ONKEY LISTENER
setReportComboOnkeypress : function () {


	// OVERRIDE LISTENER TO HANDLE KEYBOARD EVENTS
	var reportObject = this;
	dojo.connect(this.reportCombo.comboNode, "onkeypress", function(event){

		var key = event.charOrCode;

		if ( key == 13 )
		{
			reportObject.reportCombo._hideResultList();
			reportObject.newReport();
			if ( reportObject.reportCombo._popupWidget != null )
			{
				reportObject.reportCombo._showResultList();
			}
		}
	});
},



// CREATE A NEW REPORT ON TRIGGER this.reportCombo._onKeyPress ENTER
newReport : function () {

	if ( this.doingNew == true )
	{
		return;
	}

	// SET this.doingNew
	this.doingNew = true;

	// SANITY CHECK	
	if ( ! Agua.workflows )
	{
		return;
	}

	// GET NEW REPORT NAME
	var reportName = this.reportCombo.getValue();
	var projectName = this.projectCombo.getValue();
	var workflowName = this.workflowCombo.getValue();

	// CHECK IF NAME ALREADY EXISTS
	if ( Agua.isReport(projectName, reportName) )
	{
		return;
	}

	// HIDE REPORT COMBO DROPDOWN
	//this.reportCombo._isShowingNow = true;
	if ( this.reportCombo._popupWidget )
	{
		this.reportCombo._popupWidget.previousButton.style.display = "none";
		this.reportCombo._popupWidget.nextButton.style.display = "none";
		this.reportCombo._hideResultList();
	}

	// ADD REPORT TO THIS PROJECT
	Agua.addReport(
	{
		project: projectName,
		workflow: workflowName,
		report: reportName
	});

	// SHOW REPORT COMBO DROPDOWN
	this.reportCombo._showResultList();

	// UNSET this.doingNew
	this.doingNew = false;

	// RESET THE PROJECT COMBO
	this.setReportCombo(projectName, reportName);

	// SEND TO SERVER
	this.addReport(projectName, workflowName, reportName);
},


// ADD REPORT TO DATABASE ON SERVER
addReport : function (projectName, workflowName, reportName) {

	// GENERATE QUERY
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.project = projectName;
	query.workflow = workflowName;
	query.name = reportName;
	query.mode = "addReport";
	query.module = "Report::SNP";

	// POST DATA TO SERVER
	var template = this;
	dojo.xhrPost(
		{
			//url: this.url,
			url : "plugins/report/templates/SNP.json",
			contentType: "text",
			handleAs: "json",
			postData: dojo.toJson(query),
			//handleAs: "json-comment-filtered",
			timeout: 2000,
			load: function(response, ioArgs) {
				template.loadReturnValues(response, ioArgs);
			},
			error: function(response, ioArgs) {
				//return response;
			}
		}
	);
},



// DELETE A PROJECT AFTER ONCLICK 'DELETE' BUTTON BESIDE project COMBOBOX
deleteReport : function (event) {

	// SET this.doingDelete OR EXIT IF BUSY
	if ( this.doingDelete == true )
	{
		return;
	}
	this.doingDelete = true;

	if ( ! Agua.workflows )
	{
		return;
	}

	// GET DELETED PROJECT NAME OR QUIT IF EMPTY
	var projectName = this.projectCombo.getValue();
	var workflowName = this.workflowCombo.getValue();
	if ( projectName == null || ! projectName )
	{
		this.doingDelete = false;
		return;
	}

	// GET DELETED report NAME OR QUIT IF EMPTY
	var reportName = this.reportCombo.getValue();
	if ( reportName == null || ! reportName )
	{
		this.doingDelete = false;
		return;
	}

	// SET ARGS FOR CONFIRM DELETE
	var args = new Object;
	args.projectName = projectName;
	args.workflowName = workflowName;
	args.reportName = reportName;


	// DO CONFIRM DELETE
	this.confirmDelete(args);

	// UNSET this.doingDelete
	this.doingDelete = false;
},



confirmDelete : function (args) {

	var projectName = args.projectName ;
	var workflowName = args.workflowName;
	var reportName = args.reportName;

	if ( projectName == null )
	{
		return;
	}

	var reportObject = this;

	// GENERATE THE CONFIRMATION DIALOGUE
	function raiseQueryDialog(title, question, callbackFn) {

		var node = document.createElement('div');
		document.body.appendChild(node);

		var confirmDialog = new dijit.Dialog({ id: 'queryDialog', title: title }, node);
		confirmDialog.attr('class', 'deleteConfirmationDialogue');
		dojo.addClass(confirmDialog.titleBar, 'deleteConfirmationDialogueTitle');

		// When either button is pressed, kill the dialog and call the callbackFn.
		var commonCallback = function(mouseEvent)
		{
			confirmDialog.hide();
			confirmDialog.destroyRecursive();
			if (mouseEvent.explicitOriginalTarget.id == 'yesButton')
			{

				// IF THE report IS DEFINED, DELETE IT AND
				// UPDATE THE report COMBO BOX
				if ( reportName != null )
				{
					Agua.removeReport(
					{
						project: projectName,
						workflow: workflowName,
						report: reportName
					});
					reportObject.setReportCombo(projectName, workflowName);
				}
			}
			else
			{
				return false;
			}
		};

		// CREATE CONTENT AND BUTTONS
		var questionDiv = document.createElement('div');
		questionDiv.innerHTML = question;
		var yesButton = new dijit.form.Button(
			{ label: 'Yes', id: 'yesButton', onClick: commonCallback }
		);
		yesButton.attr('class', 'deleteConfirmationButton');
		dojo.addClass(yesButton, 'deleteYesButton');

		var noButton = new dijit.form.Button(
			{ label: 'No', id: 'noButton', onClick: commonCallback }
		);
		noButton.attr('class', 'deleteConfirmationButton');
		dojo.addClass(noButton, 'deleteNoButton');

		// ADD CONTENT AND BUTTONS
		dojo.addClass(yesButton, 'deleteYesButton');
		yesButton.attr('class', 'deleteConfirmationButton');

		confirmDialog.containerNode.appendChild(questionDiv);
		confirmDialog.containerNode.appendChild(yesButton.domNode);
		confirmDialog.containerNode.appendChild(noButton.domNode);

		confirmDialog.show();
	}

	var projectReport = "project '" + projectName + "'";
	if ( reportName != null )
	{
		projectReport = "report '" + reportName + "' in " + projectReport;
	}

	raiseQueryDialog("Confirm delete", "Do you really want to delete " + projectReport + "?");
},





// DO THIS LATER: TO LOAD REPORT JSON
// EXAMPLE OF xhrPost
// http://www.dojoforum.com/2007/10/11/dojo-example-xhrget-and-xhrpost
loadReport : function() {

	// GENERATE QUERY
	var query = new Object;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.project = this.projectCombo.getValue();
	query.workflow = this.workflowCombo.getValue();
	query.report = this.reportCombo.getValue();
	query.mode = "loadReport";
	query.module = "Report::SNP";

	// POST DATA TO SERVER
	var template = this;
	dojo.xhrPost(
		{
			//url: this.url,
			url : "plugins/report/templates/SNP.json",
			contentType: "text",
			handleAs: "json",
			postData: dojo.toJson(query),
			//handleAs: "json-comment-filtered",
			timeout: 2000,
			load: function(response, ioArgs) {
				template.loadReturnValues(response, ioArgs);
			},
			error: function(response, ioArgs) {
				//return response;
			}
		}
	);
},


// OPEN A WINDOW IN THE VIEW TAB WITHIN THE GIVEN CHROMOSOMAL REGION
openView : function (name, chromosome, chromosomeStart, chromosomeStop ) {



},





// OPEN FILE MANAGER TO ALLOW SELECTION OF FILE AS ARGUMENT VALUE	
openFileBrowser : function(node) {

	var valueNode = this.fileInput;

	var callArguments = new Array;
	callArguments.node = valueNode;

	// SET CALLBACK
	var selectCallback = dojo.hitch(this, function(file, location)
		{
			this.fileInput.value = file;

			//var project;
			//var workflow;
			//var filename;
			//if ( file.match( /^([^\/]+)\/([^\/]+)\/(.+)$/ ) )
			//{
			//	this.projectCombo.setValue(file.match( /^([^\/]+)\/([^\/]+)\/(.+)$/ )[1]);
			//	this.workflowCombo.setValue(file.match( /^([^\/]+)\/([^\/]+)\/(.+)$/ )[2]);
			//	this.fileInput.value = file.match( /^([^\/]+)\/([^\/]+)\/(.+)$/ )[3];
			//}
		}
	);

	// DESTROY ANY EXISTING this.fileBrowser
	if ( this.fileBrowser )
	{
		this.fileBrowser.destroy();
	}

	var fileManager = new plugins.workflow.FileManager(
	{
		paneId: "projects" + this.paneId,
		tabsNodeId: "tabs",
		project : this.project,
		workflow : this.workflow,
		selectCallback: selectCallback,
		addCallback: selectCallback,
		projects: Agua.workflows,
		workflowObject: this
	} );

	fileManager.node = node;

	this.fileBrowser = fileManager;
},




// SET EDITOR WIDGET AND ITS DEFAULT VALUE
setEditor : function () {
	this.editor = new plugins.report.Editor(
		{
			extraPlugins: "['|', 'formatBlock']",
			height :  "100px",
			minHeight: "1em",
			maxHeight: "1em"
		},
		this.editor
	);


	// SET INITIAL VALUE IF EMPTY
	var value = this.getValueType(this.editor, 'editor');


	if ( value == null )
	{
		this.setValueType(this.editor, 'editor', "A <i>rich</i> <b>text</b> summary of your report");
	}
},





// SEND ALL VALUES OF ALL ELEMENTS IN THE TEMPLATE TO THE SERVER
// editor.getValue() GIVES THE HTML MARKUP 
saveReport : function() {

	// STORE JSON QUERY
	var json = this.getData();

	// ADD mode AND class
	var query = new Object;
	query.mode = "saveReport";
	query["class"] = "Report::SNP";   

	// ADD USER AUTHENTICATION INFO AND REPORT INFO
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.project = this.projectCombo.getValue();
	query.workflow = this.workflowCombo.getValue();
	query.report = this.reportCombo.getValue();
	query.json = json;


	// POST DATA TO SERVER
	dojo.xhrPut(
		{
			url: Agua.cgiUrl + "report.cgi",
			contentType: "text",
			//handleAs: "json",
			//handleAs: "json-comment-filtered",
			putData: dojo.toJson(query),
			timeout: 20000,
			load: function(response, ioArgs) {
				return response;
			},
			error: function(response, ioArgs) {
				return response;
			}
		}
	);

},


// GET VALUES OF ALL ELEMENTS IN THE TEMPLATE
getData : function () {

	var data = new Object;
	for ( var name in this.elementObjects )
	{
		if ( ! this.elementObjects[name].valueFunction )
		{
			data[name] = '';
		}
		else
		{
			var value = this.elementObjects[name].valueFunction();
			if ( ! value )
			{
				value = false;
			}
			data[name] = value;
		}
	}

	return data;
},




// SET FILENAME IF PRESENT
setFilename : function () {

	if ( this.filename )
		this.fileInput.value= this.filename;
},


// ADD CSS CLASS INFORMATION TO COMPLETE plugins.core.ComboBox
setCombos : function () {
	// SET CSS
	this.projectCombo.popupClass = "workflow workflows dijitReset dijitMenu";
	this.projectCombo.wrapperClass = "workflows dijitPopup";
	this.projectCombo.itemHeight = 30;

	// SET CSS
	this.workflowCombo.popupClass = "workflow workflows dijitReset dijitMenu";
	this.workflowCombo.wrapperClass = "workflows dijitPopup";
	this.workflowCombo.itemHeight = 30;

	// SET CSS
	this.reportCombo.popupClass = "workflow workflows dijitReset dijitMenu";
	this.reportCombo.wrapperClass = "workflows dijitPopup";
	this.reportCombo.itemHeight = 30;

	// SET CSS
	this.speciesCombo.popupClass = "workflow workflows dijitReset dijitMenu";
	this.speciesCombo.wrapperClass = "workflows dijitPopup";
	this.speciesCombo.itemHeight = 30;

	// SET CSS
	this.chromosomeCombo.popupClass = "workflow workflows dijitReset dijitMenu";
	this.chromosomeCombo.wrapperClass = "workflows dijitPopup";
	this.chromosomeCombo.itemHeight = 30;

	// SET CSS
	this.senseCombo.popupClass = "workflow workflows dijitReset dijitMenu";
	this.senseCombo.wrapperClass = "workflows dijitPopup";
	this.senseCombo.itemHeight = 30;

	// SET CSS
	this.exonicCombo.popupClass = "workflow workflows dijitReset dijitMenu";
	this.exonicCombo.wrapperClass = "workflows dijitPopup";
	this.exonicCombo.itemHeight = 30;

	// SET CSS
	this.dbsnpCombo.popupClass = "workflow workflows dijitReset dijitMenu";
	this.dbsnpCombo.wrapperClass = "workflows dijitPopup";
	this.dbsnpCombo.itemHeight = 30;
},


//	POPULATE COMBOBOX AND SET SELECTED ITEM
setProjectCombo : function (project, workflow) {
	//console.dir(Agua);

	var projectNames = Agua.getProjectNames();

	// RETURN IF projects NOT DEFINED
	if ( ! projectNames )
	{
		return;
	}

	// SET PROJECT IF NOT DEFINED TO FIRST ENTRY IN projects
	if ( project != null || ! project)
	{
		project = projectNames[0];
	}

	// SET CSS
	this.projectCombo.popupClass = "report SNP dijitReset dijitMenu";
	this.projectCombo.wrapperClass = "SNP dijitPopup";
	this.projectCombo.itemHeight = 30;

	// SET this.project JUST IN CASE
	this.project = project;

	// DO data FOR store
	var data = {identifier: "name", items: []};
	for ( var i in projectNames )
	{
		data.items[i] = { name: projectNames[i]	};
	}

	// CREATE store
	var store = new dojo.data.ItemFileReadStore( {	data: data	} );

	//// GET PROJECT COMBO WIDGET
	if ( this.projectCombo == null )
	{
		return;
	}

	//console.dir(this.projectCombo);

	this.projectCombo.store = store;	

	// START UP AND SET VALUE
	this.projectCombo.startup();
	this.projectCombo.setValue(projectNames[0]);			

	// TO DO:
	// ADD ONKEYPRESS INSIDE COMBO INPUT
	// TO CREATE NEW PROJECTS

	// DOJO.CONNECT TO CHANGE THE this.workflowCombo
	var snpReport = this;
	dojo.connect(this.projectCombo, "onChange", function(event) {
		snpReport.project = event;
		snpReport.setWorkflowCombo(event);
	});

	// RESET THE WORKFLOW COMBO
	this.setWorkflowCombo(project, workflow);
},


// SET THE workflow COMBOBOX
setWorkflowCombo : function (project, workflow) {

	if ( project == null || ! project )
	{
		return;
	}

	// SET CSS
	this.workflowCombo.popupClass = "report SNP dijitReset dijitMenu";
	this.workflowCombo.wrapperClass = "SNP dijitPopup";
	this.workflowCombo.itemHeight = 30;


	var workflows = Agua.getWorkflowNamesByProject(project);

	// RETURN IF workflows NOT DEFINED
	if ( ! workflows )
	{
		return;
	}		

	// SET workflow IF NOT DEFINED TO FIRST ENTRY IN workflows
	if ( workflow == null || ! workflow)
	{
		workflow = workflows[0];
	}

	// SET this.workflow JUST IN CASE
	this.workflow = workflow;


	// DO data FOR store
	var data = {identifier: "name", items: []};
	for ( var i in workflows )
	{
		data.items[i] = { name: workflows[i]	};
	}

	// CREATE store
	var store = new dojo.data.ItemFileReadStore( { data: data } );

	// GET WORKFLOW COMBO
	if ( this.workflowCombo == null )
	{
		return;
	}

	this.workflowCombo.store = store;


	// GET USER INPUT WORKFLOW
	var snpReport = this;

	// SET CLASS
	//this.workflowCombo.attr('class', 'workflowWorkflowCombo');


	// START UP COMBO (?? NEEDED ??)
	this.workflowCombo.startup();
	this.workflowCombo.setValue(workflow);			

	// DOJO.CONNECT TO CHANGE THE DROP TARGET
	var snpReport = this;
	dojo.connect(this.workflowCombo, "onChange", function(event) {

		// SET this.workflow
		snpReport.workflow = event;

		// POPULATE APPLICATIONS LIST IN DROP TARGET AND INFO PANE DATA ON FIRST APPLICATION
		//var workflowApplications = snpReport.getWorkflowApplications();
	});


	this.setReportCombo(project, workflow);
},



// SET THE report COMBOBOX
setReportCombo : function (project, workflow, report) {

	if ( project == null || ! project )
	{
		return;
	}

	var reports = Agua.getReportsByProjectWorkflow(project, workflow);

	// RETURN IF reports NOT DEFINED
	if ( ! reports || reports.length == 0 )
	{
		Agua.addReport(
		{
			project: project,
			workflow: workflow,
			name: "Report1"
		});
		reports = Agua.getReportsByProjectWorkflow(project, workflow);
	}		
	var reportNames = this.hasharrayKeyToArray(reports, "name");
	if ( ! reportNames | reportNames == null )
	{
		return;
	}


	// SET report IF NOT DEFINED TO FIRST ENTRY IN reportNames
	if ( report == null || ! report)
	{
		report = reportNames[0];
	}

	// SET this.report JUST IN CASE
	this.report = report;


	// DO data FOR store
	var data = {identifier: "name", items: []};
	for ( var i in reportNames )
	{
		data.items[i] = { name: reportNames[i]	};
	}

	// CREATE store
	// http://docs.dojocampus.org/dojo/data/ItemFileWriteStore
	var store = new dojo.data.ItemFileReadStore( { data: data } );

	// GET WORKFLOW COMBO
	if ( this.reportCombo == null )
	{
		return;
	}

	this.reportCombo.store = store;

	// GET USER INPUT WORKFLOW
	var snpReport = this;

	// START UP COMBO (?? NEEDED ??)
	this.reportCombo.startup();
	this.reportCombo.setValue(report);			

	// DOJO.CONNECT TO CHANGE THE DROP TARGET
	var snpReport = this;
	dojo.connect(this.reportCombo, "onChange", function(event) {

		// SET this.report
		snpReport.report = event;

		// POPULATE APPLICATIONS LIST IN DROP TARGET AND INFO PANE DATA ON FIRST APPLICATION
		//var reportApplications = snpReport.getReportApplications();
	});
},


// SET DELETE REPORT AND FILE BROWSE BUTTONS
setButtons : function () {
	// DELETE REPORT BUTTON
	dojo.connect(this.deleteReportButton, "onclick", this, dojo.hitch( this, function(event)
		{
			this.deleteReport();
		}
	));

	// CONNECT ONCLICK
	dojo.connect(this.fileBrowse, "onclick", dojo.hitch(this, function(event)
	{
		this.openFileBrowse(event);
	}));
},
//	SET FILE UPLOAD BUTTON
setFileUpload : function () {

	// DESTROY IF EXISTS ALREADY
	if ( this.fileUploadObject != null )
	{
		console.dir(this.fileUploadObject);

		//this.fileUploadObject.destroyRecursive();

	}

	// SET NODE WHOSE VALUE WILL CHANGE TO THE UPLOADED FILE NAME
	var inputNode = this.fileInput;

	// SET NODE WHOSE VALUE WILL CHANGE TO THE UPLOADED FILE NAME
	var buttonNode = this.fileUpload;

	// SET THE PATH AS THE WORKFLOW FOLDER		
	var project = this.projectCombo.getValue();
	var workflow = this.workflowCombo.getValue();
	var path = project + "/" + workflow;

	// SET HIDDEN VALUES TO GENERATE HIDDEN NODES IN FORM
	var hiddenValues = {
		username: Agua.cookie('username'),
		sessionId : Agua.cookie('sessionId'),
		path: path
	};		

	// SET CALLBACK
	var oldFilepath = this.fileInput.value;

	var uploadUrl = Agua.cgiUrl + "upload.cgi";
	var callback = dojo.hitch(this, function(newValue)
		{				
			this.onUploadComplete(inputNode, oldFilepath, newValue);
		}
	);


	this.fileUploadObject = new plugins.upload.FileUpload(
		{
			path: path,
			inputNode: inputNode,
			buttonNode: buttonNode,
			url : uploadUrl,
			callback : callback,
			hiddenValues : hiddenValues
		}
	);
},


onUploadComplete : function (inputNode, oldValue, newValue) {

	if ( newValue )
	{
		var filename = newValue.match( /([^\/]+)$/ )[1];

		var project = this.projectCombo.getValue();
		var workflow = this.workflowCombo.getValue();
		this.fileInput.value = project + "/" + workflow + "/" + filename;	
	}

	this.setFileUpload();

},







setOkButton : function(node, name) {

	dojo.connect(this.okButton, "onclick", this, dojo.hitch( this, function(event)
		{
			// DO FILTER REPORT
			this.filterReport();
		}
	));

},


// species dijit.form.ComboBox
setSpeciesCombo : function () {

	// SET ONCHANGE TO UPDATE CHROMOSOME COMBOBOX IF CHANGED
	var snp = this;
	this.speciesCombo.onChange = function(e)
	{
		var species = e;
		snp.chromosomeSelect(species);
	}

	this.chromosomeSelect("Human");
},

// SET THE CORRECT CHROMOSOMES FOR THE SPECIES IN THE chromosome COMBOBOX
// http://en.wikipedia.org/wiki/List_of_number_of_chromosomes_of_various_organisms
chromosomeSelect : function (species) {

	var speciesChromosomes = new Object;
	speciesChromosomes["Human"] = [ "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY" ]; 
	speciesChromosomes["Mouse"] = [ "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY" ];

	// GET CHROMOSOMES FOR SPECIES
	var chromosomes = speciesChromosomes[species];

	// POPULATE chromosome SELECT BOX WITH CHROMOSOMES FOR THIS SPECIES

	// DO data FOR store
	var dataObject = {identifier: "name", items: []};
	for ( var i in chromosomes )
	{
		dataObject.items[i] = { name: chromosomes[i]	};
	}

	// CREATE store
	var storeObject = new dojo.data.ItemFileReadStore(	{	data: dataObject	} );

	// REFRESH ComboBox
	this.chromosomeCombo.store = storeObject;
	this.chromosomeCombo.startup();
	this.chromosomeCombo.setValue(chromosomes[0]);			
},

// CONNECT LISTENER TO UPDATE variantInput WHEN variantSlider CHANGES
setVariantSlider : function () {

	// ADD CLASS
	dojo.addClass(this.variantSlider.domNode, 'variant');

	// ONCHANGE
	dojo.connect(this.variantSlider, "onChange", dojo.hitch(this, function(e)
	{
		//dojo.byId(this.elementObjects['variantInput'].id).value = dojo.number.format(widget.attr('value')/100,{places:0,pattern:'#%'});
		this.variantInput.value = dojo.number.format(this.variantSlider.attr('value')/100,{places:0,pattern:'#%'});
	}));

	// ONMOUSEUP
	dojo.connect(this.variantSlider, "onMouseUp", dojo.hitch(this, function(e)
	{
		this.filterReport();
	}));

},

// INITIALISE SPINNER ON depthNode
setDepthSpinner : function () {
	this.depthNode.id = '';
	this.depthSpinner = new dijit.form.NumberSpinner(
		{ value: 10,
			constraints: {min: 1, max:1000, places:0},
			size : 1
		},
		this.depthNode
	);

	this.depthSpinner.domNode.setAttribute('style', "width: 50px;");
},

// SET CHECKBOXES TO DISABLE/ENABLE THEIR RELATED ELEMENTS
setCheckboxes : function () {

	var checkboxes = ["chromosomeCheckbox", "variantCheckbox", "depthCheckbox", "senseCheckbox", "exonicCheckbox", "dbsnp"];
	var targets = ["chromosomeCombo", "variantInput", "depthSpinner", "senseCombo", "exonicCombo", "dbsnpCombo"];	

	dojo.connect(this.chromosomeCheckbox, 'onClick', dojo.hitch(this, function ()
	{
		this.chromosomeCombo.attr('disabled', !this.chromosomeCombo.attr('disabled'));
	}));
	dojo.connect(this.variantCheckbox, 'onClick', dojo.hitch(this, function ()
	{
		// <input ... disabled="disabled">
		var disabled = this.variantInput.getAttribute('disabled');
		if ( disabled == null ) disabled = "disabled";
		else disabled = false;
		this.variantInput.disabled = disabled;
	}));
	dojo.connect(this.depthCheckbox, 'onClick', dojo.hitch(this, function ()
	{
		this.depthSpinner.attr('disabled', !this.depthSpinner.attr('disabled'));
	}));
	dojo.connect(this.senseCheckbox, 'onClick', dojo.hitch(this, function ()
	{
		this.senseCombo.attr('disabled', !this.senseCombo.attr('disabled'));
	}));
	dojo.connect(this.exonicCheckbox, 'onClick', dojo.hitch(this, function ()
	{
		this.exonicCombo.attr('disabled', !this.exonicCombo.attr('disabled'));
	}));
	dojo.connect(this.dbsnpCheckbox, 'onClick', dojo.hitch(this, function ()
	{
		this.dbsnpCombo.attr('disabled', !this.dbsnpCombo.attr('disabled'));
	}));

	// NB: THIS DOESN'T WORK
	////for ( var i = 0; i < checkboxes.length; i++ )
	//dojo.forEach(targets), function(target, i)
	//{
	//
	//	dojo.connect(this[checkboxes[i]], 'onclick', dojo.hitch(this, function ()
	//	//dojo.connect(this[checkbox], 'onclick', dojo.hitch(this, function ()
	//	{
	//		//this[targets[i]].attr('disabled', !this[targets[i]].attr('disabled'));
	//		this[target].attr('disabled', !this[target].attr('disabled'));
	//	}));
	//}
	//
},


// GET THE VALUE FOR AN OBJECT OF THE GIVEN TYPE	
getValueType : function (object, type) {
	if ( object == null )	return;
	if ( type == null )	return;

	var getSetter = this.valueFunction(type);
	if ( getSetter == null )
	{
		return;
	}

	object.valueFunction = getSetter;
	return object.valueFunction();
},


// SET THE VALUE FOR AN OBJECT OF A GIVEN TYPE
setValueType : function (object, type, value) {

	var getSetter = this.valueFunction(type);
	if ( getSetter == null )
	{
		return;
	}

	object.valueFunction = getSetter;
	return object.valueFunction(value);
},


// valueFunction
//
// PROVIDE THE GET/SET FUNCTION FOR AN OBJECT OF THE SPECIFIED TYPE
valueFunction : function (type) {

	// SET VALUE FUNCTION FOR saveReport LATER
	switch (type)
	{
		case "combobox" : case "editor" : case "checkbox" :

			return function(value) {
				if ( value )
				{
					return this.setValue(value);
				}
				else
				{
					var returnValue = this.getValue();
					if ( returnValue != null && returnValue != '' )
					{
						returnValue = returnValue.replace(/"/g, "'");
					}
					return returnValue;
				}
			};

		case "null" :
			return function(value) { return null; };

		case "radio" : return function(value) {
			return value ? this.checked = value : this.checked;
		};

		case "input" : return function(value) {
			return value ? this.value = value : this.value;
		};

		case "div" : return function(value) {
			return value ? this.innerHTML = value : this.innerHTML;
		};
		case "spinner" : return function(value) {
			return value ? this.attr('value', value) :
			dojo.number.format(this.attr('value'),{places:2}) ;
		};
	}

	return 0;
}


}); // plugins.report.SNP
