
dojo.provide("plugins.report.Template.SNP");

//// SHOW PRELOAD MESSAGE
//if ( Agua && Agua.showPreloadMessage )
//{
//	Agua.showPreloadMessage("Loading 'plugins.report.Template.SNP' ...");
//}


// USE THIS FOR LARGE DATA SETS
// http://localhost:8080/agua/html/dojo-1.5.0/dojox/grid/tests/test_sizing_ResizeHandle.html

// USE THIS FOR STYLING
// http://localhost:8080/agua/html/dojo-1.5.0/dojox/grid/tests/test_styling.html


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

// GRID
dojo.require("dojox.grid.DataGrid");

dojo.require("dojox.grid.DataGrid");
dojo.require("dojo.data.ItemFileWriteStore");
//dojo.require("dojo.parser");
dojo.require("dijit.Tooltip");
dojo.require("dijit.Menu");
dojo.require("dijit.ColorPalette");

dojo.require("dojo.parser");



// CUSTOM EDITOR
dojo.require("plugins.report.Editor");


//http://www.ensembl.org/biomart/martview/e61a39e4e8e306354f8f2a7a70bbc53c/e61a39e4e8e306354f8f2a7a70bbc53c
//http://localhost:8080/Bioptic0.2.5/html/dojo-1.5.0/demos/
//http://localhost:8080/Bioptic0.2.5/html/dojo-1.5.0/dojox/form/tests/test_SelectStack.html
//Real World Dojo part 5: Custom Components
//http://www.dotnetmafia.com/blogs/jamesashley/archive/2008/10/28/761.aspx


// FILE UPLOAD
dojo.require("plugins.upload.FileUp");


// FILE MANAGER HAS FILE SELECTORS
dojo.require("plugins.workflow.FileManager");
dojo.require("plugins.workflow.FileSelector");


// INHERITED CLASSES
dojo.require("plugins.report.Common");
dojo.require("plugins.report.Template");

dojo.declare(
    "plugins.report.Template.SNP",
    [ plugins.report.Template, plugins.report.Common ],
    {
        //Path to the template of this widget. 
        templatePath: dojo.moduleUrl("plugins", "report/templates/SNP.html"),

        // Calls dijit._Templated.widgetsInTemplate
        widgetsInTemplate : true,

        //title : "<b> SNP Report </b>",

        url : '',
	        id : '',
		filename: '',
		loading: null,

		// OR USE @import IN HTML TEMPLATE
        cssFiles : [ "dojo-1.5.0/dojox/grid/resources/Grid.css", "dojo-1.5.0/dojox/grid/resources/tundraGrid.css", "plugins/report/css/SNP.css" ],

		// CURRENT PROJECT AND WORKFLOW
		project : null,
		workflow : null,	  

		filters: [ "chromosome","variant","depth","sense","exonic","dbsnp" ],

        elements : [ "okButton", "deleteReportButton", "saveReportButton", "editor", "reportCombo", "workflowCombo", "projectCombo", "loading", "file", "totalResult", "species", "chromosomeCheckbox", "chromosomeCombo", "chromosomeResult", "variantCheckbox", "variant", "variantInput", "variantResult", "depthCheckbox", "depth", "depthResult", "sense", "senseCheckbox", "senseResult", "exonic", "exonicResult", "exonicCheckbox", "dbsnpCheckbox", "dbsnpResult", "outputResult", "fileBrowse", "fileUpload" ],

		fileBrowser : null,

        // Any initialization code would go here in the constructor.
		// plugins.report.Template and its superclasses dijit._Widget and
        // dijit._Templated do not have parameters in their constructors, so
        // there wouldn't be any multiple-inheritance complications
        // if you were to include some paramters here.
        constructor : function (args)
        {

			this.project = args.project;
			this.workflow = args.workflow;
			this.filename = args.filename;

        },



        //Inherited from dijit._Widget and called just before template
        //instantiation in buildRendering. This method is especially useful
        //for manipulating the template before it becomes visible.
        postMixInProperties: function()
        {
            //this.popup = new dijit.Dialog({});

        },





        //You can override this method to manipulate widget once it is
        //placed in the UI, but be warned that any child widgets contained
        //in it may not be ready yet.        
        postCreate: function()
        {

			// SET CLASS
			dojo.addClass(this.domNode, 'report');
			dojo.addClass(this.domNode, 'SNP');
        },



		// startup
		//
		//
		startup : function ()
		{

			// SET UP THE ELEMENT OBJECTS AND THEIR VALUE FUNCTIONS
			this.inherited(arguments);

			// NOW DO SOME MORE COMPLICATED LOGIC TO SET UP THE PAGE

			// SET PROJECT COMBO
			this.setProjectCombo();

			// SET FILE UPLOAD
			this.setFileUpload();

		},



		// setProjectCombo
		//
		//	INPUT: (OPTIONAL) project, workflow NAMES
		//
		//	OUTPUT:	POPULATE COMBOBOX AND SET SELECTED ITEM
		//
		setProjectCombo : function (project, workflow)
		{
			//console.dir(Agua);

			var projects = Agua.workflows;
			// CREATE THE DATA FOR A STORE
			if ( projects == null )
			{
				return;
			}

			if ( ! project )
			{
				// GET FIRST PROJECT
				var firstProject;

				var projectHash = projects[0];
				for ( var firstProject in projectHash )
				{
					project = firstProject;
					break;
				}
			}


			var projectNames = Agua.projectNames();

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
			var projectCombo = this.elementObjects['projectCombo'];
			if ( projectCombo == null )
			{
				return;
			}

			//console.dir(projectCombo);

			projectCombo.store = store;	

			// START UP AND SET VALUE
			projectCombo.startup();
			projectCombo.setValue(projectNames[0]);			

			// TO DO:
			// ADD ONKEYPRESS INSIDE COMBO INPUT
			// TO CREATE NEW PROJECTS

			// DOJO.CONNECT TO CHANGE THE workflowCombo
			var snpReport = this;
			dojo.connect(projectCombo, "onChange", function(event) {
				snpReport.project = event;
				snpReport.setWorkflowCombo(event);
			});

			// RESET THE WORKFLOW COMBO
			this.setWorkflowCombo(project, workflow);
		},




		// SET THE workflow COMBOBOX
		setWorkflowCombo : function (project, workflow)
		{

			if ( project == null || ! project )
			{
				return;
			}

			// CREATE THE DATA FOR A STORE
			var projectHash = Agua.getProjects();
			if ( projectHash == null )
			{
				return;
			}

			var workflows = Agua.workflowNames(project);

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
			var workflowCombo = this.elementObjects['workflowCombo'];
			if ( workflowCombo == null )
			{
				return;
			}

			workflowCombo.store = store;


			// GET USER INPUT WORKFLOW
			var snpReport = this;

			// SET CLASS
			//workflowCombo.attr('class', 'workflowWorkflowCombo');


			// START UP COMBO (?? NEEDED ??)
			workflowCombo.startup();
			workflowCombo.setValue(workflow);			

			// DOJO.CONNECT TO CHANGE THE DROP TARGET
			var snpReport = this;
			dojo.connect(workflowCombo, "onChange", function(event) {

				// SET this.workflow
				snpReport.workflow = event;

				// POPULATE APPLICATIONS LIST IN DROP TARGET AND INFO PANE DATA ON FIRST APPLICATION
				//var workflowApplications = snpReport.getWorkflowApplications();
			});


			this.setReportCombo(project, workflow);
		},




		reportCombo : function (widget, name)
		{


			//// DO data FOR store
			//var data = {identifier: "name", items: []};
			//for ( var i in projectNames )
			//{
			//	data.items[i] = { name: projectNames[i]	};
			//}


			// SET VALUE FUNCTION
			this.elementObjects[name].valueFunction = this.valueFunction('combobox');

			// OVERRIDE LISTENER TO HANDLE KEYBOARD EVENTS
			var reportObject = this;
			dojo.connect(this.elementObjects[name].comboNode, "onkeypress", function(event){

				var key = event.charOrCode;

				if ( key == 13 )
				{
					reportObject.elementObjects[name]._hideResultList();
					reportObject.newReport();
					if ( reportObject.elementObjects[name]._popupWidget != null )
					{
						reportObject.elementObjects[name]._showResultList();
					}
				}
			});
		},


		// CREATE A NEW REPORT ON TRIGGER this.reportCombo._onKeyPress ENTER
		newReport : function ()
		{

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
			var reportName = this.elementObjects['reportCombo'].valueFunction();
			var projectName = this.elementObjects['projectCombo'].valueFunction();
			var workflowName = this.elementObjects['workflowCombo'].valueFunction();

			// CHECK IF NAME ALREADY EXISTS
			if ( Agua.isReport(projectName, reportName) )
			{
				return;
			}

			// HIDE REPORT COMBO DROPDOWN
			//this.elementObjects['reportCombo']._isShowingNow = true;
			if ( this.elementObjects['reportCombo']._popupWidget )
			{
				this.elementObjects['reportCombo']._popupWidget.previousButton.style.display = "none";
				this.elementObjects['reportCombo']._popupWidget.nextButton.style.display = "none";
				this.elementObjects['reportCombo']._hideResultList();
			}

			// ADD REPORT TO THIS PROJECT
			Agua.addReport(
			{
				project: projectName,
				workflow: workflowName,
				reportName
			});

			// SHOW REPORT COMBO DROPDOWN
			this.elementObjects['reportCombo']._showResultList();

			// UNSET this.doingNew
			this.doingNew = false;

			// RESET THE PROJECT COMBO
			this.setReportCombo(projectName, reportName);


			// SEND TO SERVER
			this.addReport(projectName, workflowName, reportName);
		},







		// ADD REPORT TO DATABASE ON SERVER
		addReport : function (projectName, workflowName, reportName)
		{

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
		deleteReport : function (event)
		{

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
			var projectName = this.elementObjects['projectCombo'].valueFunction();
			var workflowName = this.elementObjects['workflowCombo'].valueFunction();
			if ( projectName == null || ! projectName )
			{
				this.doingDelete = false;
				return;
			}

			// GET DELETED report NAME OR QUIT IF EMPTY
			var reportCombo  = this.elementObjects['reportCombo'];
			var reportName = reportCombo.valueFunction();
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



		confirmDelete : function (args)
		{

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






		// SET THE report COMBOBOX
		setReportCombo : function (project, workflow, report)
		{

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
			var reportCombo = this.elementObjects['reportCombo'];
			if ( reportCombo == null )
			{
				return;
			}

			reportCombo.store = store;

			// GET USER INPUT WORKFLOW
			var snpReport = this;

			// START UP COMBO (?? NEEDED ??)
			reportCombo.startup();
			reportCombo.setValue(report);			

			// DOJO.CONNECT TO CHANGE THE DROP TARGET
			var snpReport = this;
			dojo.connect(reportCombo, "onChange", function(event) {

				// SET this.report
				snpReport.report = event;

				// POPULATE APPLICATIONS LIST IN DROP TARGET AND INFO PANE DATA ON FIRST APPLICATION
				//var reportApplications = snpReport.getReportApplications();
			});
		},






		// WORKFLOW COMBO WIDGET
		workflowCombo : function (widget, name)
		{

			// SET VALUE FUNCTION 
			this.elementObjects[name].valueFunction = this.valueFunction('combobox');
		},


		// PROJECT COMBO WIDGET
		projectCombo : function (widget, name)
		{

			// SET VALUE FUNCTION 
			this.elementObjects[name].valueFunction = this.valueFunction('combobox');
		},



		fileUpload : function (node, name)
		{

			// SET VALUE FUNCTION
			this.elementObjects[name].valueFunction = this.valueFunction("null");
		},

		setFileUpload : function ()
		{
			// SET NODE WHOSE VALUE WILL CHANGE TO THE UPLOADED FILE NAME
			var inputNode = this.elementObjects['file'];

			// SET NODE WHOSE VALUE WILL CHANGE TO THE UPLOADED FILE NAME
			var buttonNode = this.elementObjects['fileUpload'];

			// SET THE PATH AS THE WORKFLOW FOLDER		
			var project = this.elementObjects['projectCombo'].valueFunction();
			var workflow = this.elementObjects['workflowCombo'].valueFunction();
			var path = project + "/" + workflow;

			// SET HIDDEN VALUES TO GENERATE HIDDEN NODES IN FORM
			//
			// E.G.:
			//var userNode = document.createElement('input');
			//userNode.type = "hidden";
			//userNode.name = "username";
			//userNode.value = Agua.cookie('username');
			//userNode.id = "hiddenUsername";
			//_newForm.appendChild(userNode);
			//

			var hiddenValues = {
				username: Agua.cookie('username'),
				sessionId : Agua.cookie('sessionId'),
				path: path
			};		

			// SET CALLBACK
			var oldFilepath = this.elementObjects['file'].valueFunction();
			var uploadUrl = Agua.cgiUrl + "upload.cgi";
			var callback = dojo.hitch(this, function(newValue)
				{				
					this.onUploadComplete(inputNode, oldFilepath, newValue);
				}
			);



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

			//// CONNECT ONCLICK
			//dojo.connect(node, "onclick", dojo.hitch(this, function(event)
			//{
			//	this.openFileUpload(event);
			//}));
		},


		onUploadComplete : function (inputNode, oldValue, newValue)
		{

			if ( newValue )
			{
				var filename = newValue.match( /([^\/]+)$/ )[1];
				this.elementObjects['file'].valueFunction(filename);	
			}

		},



		//// OPEN FILE UPLOAD ELEMENT
		//openFileUpload : function(event)
		//{
		//	
		//	// SET THE PATH AS THE WORKFLOW FOLDER		
		//	var path = this.project + "/" + this.workflow;
		//
		//	// PASS HIDDEN VALUES TO GENERATE HIDDEN NODES IN FORM
		//	//
		//	// E.G.:
		//	//var userNode = document.createElement('input');
		//	//userNode.type = "hidden";
		//	//userNode.name = "username";
		//	//userNode.value = Agua.cookie('username');
		//	//userNode.id = "hiddenUsername";
		//	//_newForm.appendChild(userNode);
		//	//
		//
		//	var hiddenValues = {
		//		username: Agua.cookie('username'),
		//		sessionId : Agua.cookie('sessionId'),
		//		path: path
		//	};		
		//	
		//	// SET CALLBACK
		//	var valueNode = this.elementObjects['file'];
		//	var oldValue = this.elementObjects['file'].valueFunction();
		//	var uploadUrl = Agua.cgiUrl + "upload.cgi";
		//	var callback = dojo.hitch(this, function(newValue)
		//		{				
		//			this.onUploadComplete(valueNode, oldValue, newValue);
		//		}
		//	);
		//
		//	// GET THE 'BROWSE' BUTTON NODE
		//	var browseNode = this.elementObjects['fileBrowse'];
		//
		//	var fileUpload = new plugins.upload.FileUpload(
		//		{
		//			path: path,
		//			node: valueNode,
		//			browseNode : browseNode,
		//			url : uploadUrl,
		//			callback : callback,
		//			hiddenValues : hiddenValues
		//		}
		//	);
		//	
		//	////////// CREATE fileInput
		//	////////var fileInput = new dojox.form.FileInputBlind(
		//	////////	{
		//	////////		url : this.url,
		//	////////		label : 'Browse',
		//	////////		cancelText : 'Cancel',
		//	////////		blurDelay: 100,
		//	////////		triggerEvent: "onchange",
		//	////////		uploadMessage: "Uploading file..."
		//	////////	},
		//	////////	fileInputNode
		//	////////);
		//
		//	
		//	this.fileUpload = fileUpload;
		//},



		fileBrowse : function (node, name)
		{

			// SET VALUE FUNCTION
			this.elementObjects[name].valueFunction = this.valueFunction("null");

			// CONNECT ONCLICK
			dojo.connect(node, "onclick", dojo.hitch(this, function(event)
			{
				this.openFileBrowse(event);
			}));
		},



		// METHOD: openFileBrowse
		//
		// OPEN FILE MANAGER TO ALLOW SELECTION OF FILE AS ARGUMENT VALUE	
		// NB: PASS THE callback ON THROUGH FileManager, FileSelector AND _GroupSelectorPane
		//
		//
		//	FileManager:
		//		preamble: function()
		//		constructor : function(args)  ??? not necessary
		// 		loadProjectTab: function ()
		//
		//
		//	FileSelector
		//		constructor(args)		
		//		getPaneForItem: function(/*item*/ item, /* dijit._Contained */ parentPane, /* item[]? */ children){
		//f
		//
		//	_GroupSelectorPane
		//		constructor(args)		
		//		createMenu : function (type) ... dojo.connect(mItem1, "onClick", function()
		//
		//

		openFileBrowse : function(node)
		{

			var valueNode = this.elementObjects['file'];

			var callArguments = new Array;
			callArguments.node = valueNode;

			// SET CALLBACK
			var callback = dojo.hitch(this, function(file, location)
				{

					var project;
					var workflow;
					var filename;
					if ( file.match( /^([^\/]+)\/([^\/]+)\/(.+)$/ ) )
					{
						this.elementObjects['projectCombo'].valueFunction(file.match( /^([^\/]+)\/([^\/]+)\/(.+)$/ )[1]);
						this.elementObjects['workflowCombo'].valueFunction(file.match( /^([^\/]+)\/([^\/]+)\/(.+)$/ )[2]);
						this.elementObjects['file'].valueFunction(file.match( /^([^\/]+)\/([^\/]+)\/(.+)$/ )[3]);
					}
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
				callback: callback,
				projects: Agua.workflows,
				workflowObject: this
			} );

			fileManager.node = node;

			this.fileBrowser = fileManager;
		},




		// editor
		//
		// NOTES ABOUT THE REPORT

        editor : function (node, name)
        {

			var widget = new plugins.report.Editor(
				{
					extraPlugins: "['|', 'formatBlock']",
					height :  "100px",
					minHeight: "1em",
					maxHeight: "1em"
				},
				node
			);

			// RESET OBJECT TO WIDGET
			this.elementObjects[name] = widget;
			this.elementObjects[name].valueFunction = this.valueFunction('editor');

			// SET INITIAL VALUE IF EMPTY
			if ( ! this.elementObjects[name].valueFunction || this.elementObjects[name].valueFunction() == '' )
			{
				this.elementObjects[name].valueFunction("A <i>rich</i> <b>text</b> summary of your report");
			}
        },




		filterReport : function()
		{
            var snpReport = this;

			// SET AND DISPLAY LOADING
            this.loading = true;
            var loadingIcon = this.elementObjects["loading"];
            dojo.removeClass(loadingIcon, "static");
            dojo.addClass(loadingIcon, "loading");

			// GENERATE QUERY
			var query = new Object;
			query.username = Agua.cookie('username');
			query.sessionId = Agua.cookie('sessionId');
			query.project = Agua.cookie('project');
			query.workflow = Agua.cookie('workflow');
			query.report = Agua.cookie('report');

			// GET THE CURRENT SETTINGS OF ALL THE REPORT FORM WIDGETS
			var query = new Object;
			for ( name in this.elementObjects )
			{
				if ( this.elementObjects[name].valueFunction )
				{
					var value;
					if ( name == "file" )
					{
						value = this.elementObjects['projectCombo'].valueFunction();
						value += "/";
						value += this.elementObjects['workflowCombo'].valueFunction();
						value += "/";
						value += this.elementObjects[name].valueFunction();
					}
					else
					{
						value = this.elementObjects[name].valueFunction();
					}

					if ( ! value )
					{
						value = false;
					}
					query[name] = value;
				}
			}

			// ADD mode AND class
			query.mode = "filterReport";
			query.module = "Report::SNP";
			query.username = Agua.cookie('username');
			query.sessionId = Agua.cookie('sessionId');
			query.project = Agua.cookie('project');
			query.workflow = Agua.cookie('workflow');
			query.report = Agua.cookie('report');

			// GET THE FILTER VALUES
			var filters = this.filters;
            var filtersObject = new Object;
			for ( var index in filters )
			{
				var name = filters[index] + "Checkbox";
				var value;
				if ( this.elementObjects[name].valueFunction )
				{
					var value = this.elementObjects[name].valueFunction();
					if ( value == "on" )	{	value = true;	}
					if ( value == "off" )	{	value = false;	}
					filtersObject[filters[index]] = value;
				}
			}
			query.filters = filtersObject;

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
                        snpReport.loadReturnValues(response, ioArgs);
                    },
					error: function(response, ioArgs) {
						return response;
					}
				}
			);
		},


		// loadReturnValues
		//
		// LOAD RETURNED VALUES FROM SERVER INTO WIDGETS/NODES ON PAGE
		//
        loadReturnValues : function (filters)
        {

			// SET UP FILTERS RESULTS IN ORDER
			var filterNames = this.filters;
			filterNames.unshift("total");
			filterNames.push("output");
			for ( var i = 0; i < filterNames.length; i++ )
			{
				filterNames[i] = filterNames[i] + "Result";
			}

            for ( var i = 0; i < filterNames.length; i++)
			{
				var filter  = filterNames[i];

				if ( filter == "outputResult" )
				{
					var data = {"identifier":"id","label":"id","items":[]};

					// CONVERT ARRAYS INTO HASHES
					var fieldNames = ['name','chromosome','ccdsstart','ccdsstop','referencenucleotide','variantnucleotide','depth','variantfrequency','chromosomestart','chromosomestop','sense','referencecodon','variantcodon','referenceaa','variantaa','strand','snp','score','dbsnpstrand'];
					var twoDarray = dojo.fromJson(filters[filter]);
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
	                this.elementObjects[filter].valueFunction(data);
				}
				else
				{
	                this.elementObjects[filter].valueFunction(filters[filter]);
				}
			}

            // DISPLAY FINISHED LOADING
            this.loading = false;
            var loadingIcon = this.elementObjects["loading"];
            dojo.removeClass(loadingIcon, "loading");
            dojo.addClass(loadingIcon, "static");
        },


		okButton : function(node, name)
		{

			// SET VALUE FUNCTION 
			this.elementObjects[name].valueFunction = this.valueFunction('null');

			dojo.connect(node, "onclick", this, dojo.hitch( this, function(event)
				{
					// DO FILTER REPORT
					this.filterReport();
				}
			));

		},





		saveReportButton : function(node, name)
		{

			// SET VALUE FUNCTION 
			this.elementObjects[name].valueFunction = this.valueFunction('null');

			dojo.connect(node, "onclick", this, dojo.hitch( this, function(event)
				{
					// SAVE REPORT
					this.saveReport();
				}
			));

		},


		deleteReportButton : function(node, name)
		{

			// SET VALUE FUNCTION 
			this.elementObjects[name].valueFunction = this.valueFunction('null');

			dojo.connect(node, "onclick", this, dojo.hitch( this, function(event)
				{
					// delete REPORT
					this.deleteReport();
				}
			));

		},





		loadReportButton : function(widget, name)
		{

			var snp = this;
			widget.onClick = function(e)
			{
				snp.loadReport();
			}

			// SET VALUE FUNCTION FOR loadReport LATER
			this.elementObjects[name].valueFunction = this.valueFunction('null');
		},



		// EXAMPLE OF xhrPost
        // http://www.dojoforum.com/2007/10/11/dojo-example-xhrget-and-xhrpost

		loadReport : function()
		{

			// GENERATE QUERY
			var query = new Object;
			query.username = Agua.cookie('username');
			query.sessionId = Agua.cookie('sessionId');
			query.project = Agua.cookie('project');
			query.workflow = Agua.cookie('workflow');
			query.report = Agua.cookie('Report1');
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




        loading : function (node, name){

			// ADD 'HIDE' CLASS
			//dojo.addClass(node, 'hide');

			// VALUE FUNCTION 
			this.elementObjects[name].valueFunction = this.valueFunction('null');
        },



        outputResult : function(widget, name)
		{

            var snpReport = this;

			// SET VALUE FUNCTION
			this.elementObjects[name].valueFunction = function(value)
			{

                if ( value )
                {
                    snpReport.loadGrid(value);
                }                
            }			
		},




        loadGrid : function(data)
        {

            // REMOVE EXISTING GRID NODE FROM GRID CONTAINER
            var gridContainer = dojo.byId(this.elementObjects["outputResult"].id);
            while(gridContainer.firstChild)
            {
               gridContainer.removeChild(gridContainer.firstChild);
            }

            // SET GRID ID AND REMOVE ANY EXISTING GRID
            var gridId = this.gridId;
            if ( dijit.byId(gridId) )
            {
                dijit.byId(gridId).destroy();
            }

			// SEE EXAMPLE:    
            //http://localhost:8080/Bioptic0.2.5/html/tests/report-snp3.html

            // SET STORE USING DATA
            var gridStore = new dojo.data.ItemFileWriteStore({data: data});


            var layout = [
                {name: 'Entry', field: 'id', width: "30px"},
                {name: 'Name', field: 'name', width: "100px"},
                {name: 'Chr', field: 'chromosome', width: "40px"},
                {name: 'CCDS start', field: 'ccdsstart', width: "40px"},
                {name: 'CCDS stop', field: 'ccdsstop', width: "40px"},
                {name: 'Ref nt', field: 'referencenucleotide', width: "20px"},
                {name: 'Var nt', field: 'variantnucleotide', width: "20px"},
                {name: 'Depth', field: 'depth', width: "30px"},
                {name: 'Var freq', field: 'variantfrequency', width: "30px"},
                {name: 'Chr start', field: 'chromosomestart', width: "40px"},
                {name: 'Chr stop', field: 'chromosomestop', width: "40px"},
                {name: 'Sense', field: 'sense', width: "40px"},
                {name: 'Ref codon', field: 'referencecodon', width: "35px"},
                {name: 'Var codon', field: 'variantcodon', width: "35px"},
                {name: 'Ref aa', field: 'referenceaa', width: "50px"},
                {name: 'Var aa', field: 'variantaa', width: "50px"},
                {name: 'Strand', field: 'strand', width: "25px"},
                {name: 'SNP', field: 'snp', width: "40px"},
                {name: 'Score', field: 'score', width: "40px"},
                {name: 'dbSNP strand', field: 'dbsnpstrand', width: "40px"}
            ];

            var gridId = dojo.dnd.getUniqueId();

            grid = new dojox.grid.DataGrid({
                id: gridId,
				autoHeight: true,
                store: gridStore,
                structure: layout
            }, document.createElement('div'));

            // IMPORTANT: CLASS IS REQUIRED FOR PROPER DISPLAY
            dojo.addClass(grid.viewsHeaderNode, 'viewsHeaderNode');
            dojo.addClass(grid.domNode, 'reportSNPGrid');
			dojo.addClass(gridContainer, 'work-pane-split');

//
//            var dojoxGridViewNode = grid.viewsNode.firstChild;
//			//dojoxGridViewNode children
//			//0
//			//	input.dojoxGridHiddenFocus on
//			//1
//			//	input.dojoxGridHiddenFocus on
//			//2
//			//	div.dojoxGridScrollbox
//            var dojoxGridScrollbox = dojoxGridViewNode.childNodes[2];


			//dojo.addClass(grid.headerContentNode, 'headerContentNode');



            // APPEND GRID TO GRID CONTAINER AND STARTUP
            gridContainer.appendChild(grid.domNode);

            grid.startup();

            // RESPONDS TO CLICKS ON THE GRID -- BUT CAN'T FIND GRID...
            this.grid = grid;

            var snpReport = this;
            dojo.connect(grid,'onClick', function getGridRow (n)
                {
                    var grid = snpReport.grid;

                    var row = grid.selection.getSelected()[0];

                    var datarow = grid.getItem(n.rowIndex);


                    var name = grid.store.getValue(datarow,"name");
                    var chromosome = grid.store.getValue(datarow,"chromosome");
                    var chromosomeStart = grid.store.getValue(datarow,"chromosomestart");
                    var chromosomeStop = grid.store.getValue(datarow,"chromosomestop");


                    //var t = grid.store.getValue(row,"Title");
                    //console.debug("jjjj jjjj f: " +f+" "+t);

                    return;
                }
            );



//            // ADD CONTEXT MENU
            var gridMenu = this.createGridMenu(grid);


////          // ADD TOOLTIP
////            var gridTooltipEnabled = true;
////			var showTooltip = function(e) {
////
////
////                if( gridTooltipEnabled )
////                {
////                    var msg = "This is cell " + e.rowIndex + ", " + e.cellIndex;
////                    dijit.showTooltip(msg, e.cellNode);
////                }
////                
////            };
////
////            var hideTooltip = function(e) {
////                dijit.hideTooltip(e.cellNode);
////                // FIXME: make sure that pesky tooltip doesn't reappear!
////                // would be nice if there were a way to hide tooltip without regard to aroundNode.
////                dijit._masterTT._onDeck=null;
////            };
////			
////            // cell tooltip
////            dojo.connect(grid, "onClick", showTooltip);
////            dojo.connect(grid, "onClick", hideTooltip);
////            //dojo.connect(grid, "onCellMouseOver", showTooltip);
////            //dojo.connect(grid, "onCellMouseOut", hideTooltip);
////
////            // header cell tooltip
////            dojo.connect(grid, "onHeaderCellMouseOver", showTooltip);
////            dojo.connect(grid, "onHeaderCellMouseOut", hideTooltip);


            // grid menu
            //window["gridMenu"] = dijit.byId("gridMenu");



// COMMENTED OUT FOR DEBUG			
			//gridMenu.bindDomNode(grid.domNode);





            //// prevent grid methods from killing the context menu event by implementing our own handler
            //grid.onCellContextMenu = function(e) {
            //    cellNode = e.cellNode;
            //};
            //grid.onHeaderContextMenu = function(e) {
            //    cellNode = e.cellNode;
            //};


        },   // loadGrid




        openView : function (name, chromosome, chromosomeStart, chromosomeStop )
        {





        },



        // ADD PROGRAMMATIC CONTEXT MENU
        createGridMenu : function (grid)
        {
            var menuId = dojo.dnd.getUniqueId();
            var dynamicMenu = new dijit.Menu( { id: menuId } );

            // ADD MENU TITLE
            dynamicMenu.addChild(new dijit.MenuItem( { label:"Options", disabled:false} ));
            dynamicMenu.addChild(new dijit.MenuSeparator());

            //// ONE OF FOUR WAYS TO DO MENU CALLBACK WITH ACCESS TO THE MENU ITEM AND THE CURRENT TARGET 	
            // 4. dojo.connect CALL
            //	REQUIRES:
            //		ADDED menu.currentTarget SLOT TO dijit.menu
            var mItem1 = new dijit.MenuItem(
                {
                    label: "View",
                    disabled: false
                }
            );
            dynamicMenu.addChild(mItem1);


            var snpReport =this;
            dojo.connect(mItem1, "onClick", function(e)
                {


                    var row = grid.selection.getSelected()[0];

                    //var datarow = grid.getItem(row.rowIndex);


                    var currentTarget = e.currentTarget;
                    var rowIndex = e.target.rowIndex;

                    //var row = grid.selection.getSelected()[0];


            return;

                    snpReport.viewRow(rowIndex, grid);

                }
            );

            // SEPARATOR
            dynamicMenu.addChild(new dijit.MenuSeparator());

            //	ADD run MENU ITEM
            var mItem2 = new dijit.MenuItem(
                {
                    label: "Edit",
                    disabled: false
                }
            );
            dynamicMenu.addChild(mItem2);	

            dojo.connect(mItem2, "onClick", function()
                {

                    var currentTarget = dynamicMenu.currentTarget; 
                    var adminList = currentTarget.parentNode;
                }
            );

            return dynamicMenu;

        },   // createGridMenu


        // DISPLAY THE GRID ITEM IN THE 'View' TAB
        viewRow : function (rowIndex, grid)
        {
            //var value = row[0];

            var row = grid.selection.getSelected()[0];
            //var datarow = grid.getItem(rowIndex);
            var name = grid.store.getValue(datarow,"name");
            var chromosome = grid.store.getValue(datarow,"chromosome");
            var chromosomeStart = grid.store.getValue(datarow,"chromosomestart");
            var chromosomeStop = grid.store.getValue(datarow,"chromosomestart");

            //var t = grid.store.getValue(row,"Title");
            //console.debug(f+" "+t);

            //var store = grid.store;

            //var row = store.fetch(
            //{
            //    query : { id : rowIndex }, 
            //    onComplete: function(items, request)
            //    {
            //        
            //        return 1;
            //    }
            //});

            //var row = grid.getItem(rowIndex);
            //
            //var value = row[0];
            //
            //var row = grid.selection.getSelected()[rowIndex];
            //console.debug("row:" + row);

            //var value = grid.model.data[row];
            //console.debug("value: " + value);

            //var erName = value.ER_Name;
            //console.debug("erName : " + erName);


        },



		file : function (node, name)
		{


			// SET VALUE FUNCTION
			this.elementObjects[name].valueFunction = this.valueFunction('textinput');

			// SET FILENAME IF PRESENT
			if ( this.filename )
			{
				this.elementObjects[name].valueFunction(this.filename);
			}
		},




		totalResult : function (node, name)
		{

			// SET VALUE FUNCTION 
			this.elementObjects[name].valueFunction = this.valueFunction('div');
		},



		// species dijit.form.ComboBox
		//
		//

		species : function (widget, name)
		{

			// SET ONCHANGE TO UPDATE CHROMOSOME COMBOBOX IF CHANGED
			var snp = this;
			widget.onChange = function(e)
			{
				var species = e;
				snp.chromosomeSelect(species);
			}

			// SET VALUE FUNCTION
			this.elementObjects[name].valueFunction = this.valueFunction('combobox');
		},



		chromosomeCheckbox : function (widget, name)
		{

			dojo.connect(widget, 'onClick', dojo.hitch(this, function ()
			{
				this.elementObjects['chromosomeCombo'].attr('disabled', !this.elementObjects['chromosomeCombo'].attr('disabled'));
			}));

			// SET VALUE FUNCTION
			this.elementObjects[name].valueFunction = this.valueFunction('checkbox');
		},


		// chromosomeCombo
		//
		// chromosome dijit.form.ComboBox		
        // http://en.wikipedia.org/wiki/List_of_number_of_chromosomes_of_various_organisms
		//

		chromosomeCombo : function (widget, name)
		{

			// SET VALUE FUNCTION
			this.elementObjects[name].valueFunction = this.valueFunction('combobox');

			dojo.require("dojo.data.ItemFileReadStore");			
			this.chromosomeSelect('Human');
		},



		// SET THE CORRECT CHROMOSOMES FOR THE SPECIES IN THE chromosome COMBOBOX
		// http://en.wikipedia.org/wiki/List_of_number_of_chromosomes_of_various_organisms
		chromosomeSelect : function (species)
		{

			var speciesChromosomes = new Object;
			speciesChromosomes["Human"] = [ "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr67", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY" ]; 
			speciesChromosomes["Mouse"] = [ "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr67", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY" ];

			// GET CHROMOSOMES FOR SPECIES
			var chromosomes = speciesChromosomes[species];

			// POPULATE chromosome SELECT BOX WITH CHROMOSOMES FOR THIS SPECIES
			var chromosomeCombobox = this.elementObjects["chromosomeCombo"];

			// DO data FOR store
			var dataObject = {identifier: "name", items: []};
			for ( var i in chromosomes )
			{
				dataObject.items[i] = { name: chromosomes[i]	};
			}

			// CREATE store
			var storeObject = new dojo.data.ItemFileReadStore(	{	data: dataObject	} );

			// REFRESH ComboBox
			chromosomeCombobox.store = storeObject;
			chromosomeCombobox.startup();
			chromosomeCombobox.setValue(chromosomes[0]);			
		},


		chromosomeResult : function (widget, name)
		{
			// SET VALUE FUNCTION 
			this.elementObjects[name].valueFunction = this.valueFunction('div');
		},


        variantCheckbox : function (widget, name)
        {

            var snp = this;
			widget.onChange = function()
			{
				var variant = snp.elementObjects["variant"];
				variant.attr('disabled', !variant.attr('disabled'));
			}

			// SET VALUE FUNCTION
			this.elementObjects[name].valueFunction = this.valueFunction('checkbox');
		},


        variant : function (widget, name)
        {

			// ADD CLASS
			dojo.addClass(widget.domNode, 'variant');

			// ONCHANGE
			dojo.connect(widget, "onChange", dojo.hitch(this, function(e)
			{
				dojo.byId(this.elementObjects['variantInput'].id).value = dojo.number.format(widget.attr('value')/100,{places:0,pattern:'#%'});
			}));

			// ONMOUSEUP
			dojo.connect( widget, "onMouseUp", dojo.hitch(this, function(e)
			{
				this.filterReport();
			}));

			// SET VALUE FUNCTION
			this.elementObjects[name].valueFunction = function(value)
			{
				return value ? this.attr('value', value) :
                dojo.number.format(this.attr('value'), { places: 0 });
			}
        },


		variantInput : function (node, name)
		{

			dojo.addClass(node, 'variantInput');

			// SET VALUE FUNCTION
			this.elementObjects[name].valueFunction = this.valueFunction('textinput');
		},


		variantResult : function (widget, name)
		{
			// SET VALUE FUNCTION
			this.elementObjects[name].valueFunction = this.valueFunction('div');
		},


		depthCheckbox : function (widget, name)
		{
			var snpReport = this;
			widget.onClick = function (e)
			{
				var depth = snpReport.getWidget('depth');
				if ( depth )
				{
					depth.attr('disabled', !depth.attr('disabled'));
				}
			}

			// SET VALUE FUNCTION
			this.elementObjects[name].valueFunction = this.valueFunction('checkbox');
		},



		depth : function (node, name)
		{
			node.id = '';
			var widget = new dijit.form.NumberSpinner(
				{ value: 10,
					constraints: {min: 1, max:1000, places:0},
					size : 1
				},
				node
			);

			this.elementObjects[name] = widget;
			dojo.addClass(widget.domNode, 'depth');

			// SET VALUE FUNCTION
			this.elementObjects[name].valueFunction = function(value)
			{
				return value ? this.attr('value', value) :
					dojo.number.format(this.attr('value'),{places:2}) ;
			}
		},


		depthResult : function (widget, name)
		{
			// SET VALUE FUNCTION 
			this.elementObjects[name].valueFunction = this.valueFunction('div');
		},


		senseCheckbox : function (widget, name)
		{
			var snpReport = this;
			widget.onClick = function (e)
			{
				var sense = snpReport.elementObjects['sense'];
				if ( sense )
				{
					sense.attr('disabled', !sense.attr('disabled'));
				}
			}

			// SET VALUE FUNCTION 
			this.elementObjects[name].valueFunction = this.valueFunction('radio');
		},


		sense : function (widget, name)
		{
			// SET VALUE FUNCTION 
			this.elementObjects[name].valueFunction = this.valueFunction('combobox');
		},


		senseResult : function (widget, name)
		{
			// SET VALUE FUNCTION 
			this.elementObjects[name].valueFunction = this.valueFunction('div');
		},


		exonicCheckbox : function (widget, name)
		{


			var snpReport = this;
			widget.onClick = function (e)
			{
				var exonic = snpReport.elementObjects['exonic'];
				if ( exonic )
				{
					exonic.attr('disabled', !exonic.attr('disabled'));
				}
			}

			// SET VALUE FUNCTION 
			this.elementObjects[name].valueFunction = this.valueFunction('radio');
		},


		exonic : function (widget, name)
		{
			// SET VALUE FUNCTION 
			this.elementObjects[name].valueFunction = this.valueFunction('combobox');
		},


		exonicResult : function (widget, name)
		{
			// SET VALUE FUNCTION 
			this.elementObjects[name].valueFunction = this.valueFunction('div');
		},


		dbsnpCheckbox : function (widget, name)
		{
			// SET VALUE FUNCTION 
			this.elementObjects[name].valueFunction = this.valueFunction('radio');
		},


		dbsnpResult : function (widget, name)
		{
			// SET VALUE FUNCTION 
			this.elementObjects[name].valueFunction = this.valueFunction('div');
		},



		// output
		//
		// GRID REPORT GOES HERE

        output : function (widget, name)
        {
			var snpReport = this;

			// SET VALUE FUNCTION
			if ( ! this.elementObjects[name] )
			{
			}



			this.elementObjects[name].valueFunction = function(widgetObject, value)
			{
				if ( widgetObject.widget )
				{
					return value ? widgetObject.widget.setValue(value) : widgetObject.widget.getValue();
				}
				else if ( dojo.byId(widgetObject.id ) && ( value == null || ! value ) )
				{
					return dojo.byId(widgetObject.id).value; 
				}
			}
        }
	}

); // plugins.report.Template.SNP



        // NOT NEEDED AS USED dojo.connect DIRECTLY TO GRID FOR ONCLICK

        //getGridRow : function (n)
        //{
        //    var grid = this.grid;
        //    
        //    var row = grid.selection.getSelected()[0];
        //    
        //    var datarow = grid.getItem(n.rowIndex);
        //    
        //    var f = grid.store.getValue(datarow,"name");
        //    var t = grid.store.getValue(row,"Title");
        //    console.debug("jjjj jjjj f: " +f+" "+t);
        //
        //    return;
        //},





}
