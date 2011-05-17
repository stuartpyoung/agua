dojo.provide("plugins.workflow.Workflow");

/*	CLASS SUMMARY

	THIS IS THE MAIN CLASS IN THE WORKFLOW PLUGIN WHICH INSTANTIATES ALL THE OTHERS

	1. SET LEFT PANE
		Apps.js
		SharedApps.js

	2. MIDDLE PANE
		Stages.js --> Target.js/StageRow.js


	3. RIGHT PANE
		SET STAGES, HISTORY AND SHARED IN 



	WIDGET HIERARCHY

				1_to_1	   1_to_1	 
		Workflow --> Stages --> Target 
		 |              |
		 |              |1_to_many
		 |              |   
		 | 1_to_1       --> StageRow
		 |                     |
		 |                     | 1_to_1
		 |                     |
		 -----------------> Parameters --> ParameterRow
                                    1_to_many

    CORE LIST

    core.workflow       =   plugins.workflow.Workflow
    core.parameters  =   plugins.workflow.Parameters
    core.stages         =   plugins.workflow.Stages
    core.target         =   plugins.workflow.Target
    core.fileManager    =   plugins.workflow.FileManager
    core.apps           =   plugins.workflow.Apps
    core.sharedApps     =   plugins.workflow.SharedApps
    core.sharedProjects =   plugins.workflow.SharedProjects


	USAGE SCENARIO 1: CREATION AND LOADING OF A NEW WORKFLOW PANE

		--> CREATE RIGHT PANE Parameters.js AS this.Parameters

		--> CREATE MIDDLE PANE Stages.js VIA ITS METHOD setDropTarget

			--> CREATE DROP TARGET Target.js

				--> OVERRIDE onDndDrop TO CONVERT DROPPED NODE

					INTO StageRow.js WITH ONCLICK loadParametersPane

					(CALLS loadParametersPane METHOD IN Workflow.js

					WHICH IN TURN CALLS load METHOD OF Parameters.js

			--> CALL loadParametersPane METHOD IN Workflow.js

				--> CALLS load METHOD OF Parameters.js

					--> CALLS checkValidParameters IN StageRow

			--> CHECK VALIDITY OF OTHER StageRows (2, 3, 4, ...)

				--> CALL checkValidParameters IN StageRow

			--> UPDATE VALIDITY OF Stage.js (RunWorkflow BUTTON)


	USAGE SCENARIO 2: USER DROPS APPLICATION INTO TARGET

		1. onDndDrop METHOD IN Target.js

			--> CONVERTS App.js INTO StageRow.js

			--> CALLS loadParametersPane IN Workflow.js 

				--> CALLS load IN Parameters.js

					--> CALLS checkValidParameters IN StageRow


	USAGE SCENARIO 3: USER UPDATES PARAMETER IN DATA PANE

		1. ParameterRow.js CHECKS VALIDITY AND PRESENCE OF FILES

			--> CALLS checkValidInputs METHOD OF Parameters.js

				GETS this.isValid FROM VALIDITY OF ALL PARAMETERS 

				--> CALLS setValid/setInvalid OF StageRow.js 

					--> CALLS updateValidity OF Stages.js

						POLL VALIDITY OF ALL StageRows

						SET RunWorkflow BUTTON IF ALL STAGES VALID

*/


// EXPANDOPANE
dojo.require("dojox.layout.ExpandoPane");

// STORE FOR COMBO BOXES
dojo.require("dojo.data.ItemFileReadStore");

// FILE UPLOAD
dojo.require("plugins.upload.FileUpload");

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

// FILE MANAGER & FILE SELECTORS
dojo.require("plugins.workflow.FileManager");
dojo.require("plugins.files._FileInfoPane");
dojo.require("plugins.workflow.FileSelector");

// WIDGETS AND TOOLS FOR EXPANDO PANE
dojo.require("dojo.data.ItemFileReadStore");
dojo.require("dijit.form.ComboBox");
dojo.require("dijit.Tree");
dojo.require("dijit.layout.AccordionContainer");
dojo.require("dijit.layout.TabContainer");
dojo.require("dijit.layout.ContentPane");
dojo.require("dijit.layout.BorderContainer");
dojo.require("dojox.layout.FloatingPane");
dojo.require("dojo.fx.easing");
dojo.require("dojox.rpc.Service");
dojo.require("dojo.io.script");
dojo.require("dijit.TitlePane");

// DnD
dojo.require("dojo.dnd.Source"); // Source & Target
dojo.require("dojo.dnd.Moveable");
dojo.require("dojo.dnd.Mover");
dojo.require("dojo.dnd.move");

// Menu
dojo.require("plugins.menu.Menu");

// TIMER
dojo.require("dojox.timing");

// TOOLTIP
dojo.require("dijit.Tooltip");

// TOOLTIP DIALOGUE
dojo.require("dijit.Dialog");
dojo.require("dijit.form.Textarea");
dojo.require("dijit.Editor");
dojo.require("dijit.form.CheckBox");
dojo.require("dijit.form.Button");

// INHERITED
dojo.require("plugins.core.Common");

// LAYOUT WIDGETS
dojo.require("dijit.layout.SplitContainer");
dojo.require("dijit.layout.ContentPane");

// PANES
dojo.require("plugins.core.BorderContainerStatic");
dojo.require("plugins.core.BorderContainer");
dojo.require("plugins.core.ExpandoPane");

// RUN STATUS
dojo.require("plugins.workflow.RunStatus");

// HISTORY
dojo.require("plugins.workflow.History");
dojo.require("plugins.workflow.HistoryPane");

// SHARED
dojo.require("plugins.workflow.SharedApps");
dojo.require("plugins.workflow.SharedProjects");

// INFO PANE (DISPLAY STAGE DATA)
dojo.require("plugins.workflow.Parameters");

// SHARED APP ROW
dojo.require("plugins.workflow.AppRow");

// APPS DRAG SOURCE
dojo.require("plugins.workflow.Apps");

// WORKLOWS DROP TARGET
dojo.require("plugins.workflow.Stages");


dojo.declare( "plugins.workflow.Workflow",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "workflow/templates/workflow.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

// CSS FILE FOR BUTTON STYLING
cssFiles : [
    dojo.moduleUrl("plugins") + "/workflow/css/workflow.css",
    dojo.moduleUrl("plugins") + "/workflow/css/history.css",
    dojo.moduleUrl("plugins") + "/workflow/css/shared.css"
],

// PARENT NODE, I.E., TABS NODE
attachWidget : null,

// PROJECT NAME AND WORKFLOW NAME IF AVAILABLE
project : null,
workflow : null,

// POLL SERVER FOR WORKFLOW STATUS
polling : false,

// INSERT TEXT BREAKS WIDTH, CORRESPONDS TO CSS WIDTH OF INPUT 'value' TABLE ELEMENT
textBreakWidth : 22,

// plugins.workflow.FileManager
fileManager : null,

// CORE WORKFLOW OBJECTS
core : new Object,

	////}}

constructor : function(args) {
// CONSTRUCTOR	
// Any initialization code would go here in the constructor.
// plugins.report.Template and its superclasses dijit._Widget and
// dijit._Templated do not have parameters in their constructors, so
// there wouldn't be any multiple-inheritance complications
// if you were to include some paramters here.

	// LOAD CSS
	this.loadCSS();

	// SET ARGS
	this.attachWidget = Agua.tabs;
	if ( args != null )
	{
		this.project = args.project;
		this.workflow = args.workflow;
	}

	// SET CORE CLASSES
	this.core.workflow = this;
},

postMixInProperties: function() {
//Inherited from dijit._Widget and called just before template
//instantiation in buildRendering. This method is especially useful
//for manipulating the template before it becomes visible.
	//this.popup = new dijit.Dialog({});

},

postCreate: function() {
//You can override this method to manipulate widget once it is
//placed in the UI, but be warned that any child widgets contained
//in it may not be ready yet.        

	// DON'T NEED THIS AS CLASS IS SET IN TEMPLATE
	// SET CLASS
	//dojo.addClass(this.domNode, 'workflow');
	//dojo.addClass(this.domNode, 'SNP');

	this.startup();
},


startup : function () {
// START UP CHILD PANES

//// START STANDBY
//dojo.require("dojox.widget.Standby");
//
//        var standby = new dojox.widget.Standby({
//            target: Agua.tabs.containerNode
//        });
//		document.body.appendChild(standby.domNode);
//
//        standby.show();
//

	// SET UP THE ELEMENT OBJECTS AND THEIR VALUE FUNCTIONS
	this.inherited(arguments);

    // ADD THIS WIDGET TO Agua.workflowWidgets
    Agua.addWidget("workflow", this);

	// ADD THE PANE TO THE TAB CONTAINER
	this.attachWidget.addChild(this.mainTab);
	this.attachWidget.selectChild(this.mainTab);

	// EXPAND LEFT PANE (WAS CLOSED SO THAT RIGHT PANE WOULD RENDER)
	this.leftPaneExpando.toggle();

	// INSTANTIATE PARAMETERS PANE (RIGHT PANE)
	this.setParameters();

	//// SET APPLICATIONS (LEFT PANE)
	//this.setApps();

	// SET SHARED APPLICATIONS (LEFT PANE)
	this.setSharedApps();

	// INSTANTIATE FILE MANAGER
	this.setFileManager();

	// SET WORKFLOWS DROP TARGET (MIDDLE PANE)
	this.setStages();

	//// SET SHARED PROJECTS
	//this.setSharedProjects();	

	//// SET HISTORY PANE AND GET WORKFLOW HISTORY
	//this.setHistory();
	//
	//// SET Stage AS SELECTED IN MIDDLE PANE
	//this.middlePane.selectChild(this.stages.mainTab);

	// INSTANTIATE RUN STATUS (RIGHT PANE)
	this.setRunStatus();

},


///////////////////////////////////////////////////////////////////////////
/////////////////		SET TABS IN PANES
setFileManager : function() {
/* OPEN FILE MANAGER TO ALLOW SELECTION OF FILE AS ARGUMENT VALUE	
	NB: PASS THE callback ON THROUGH FileManager, FileSelector AND _GroupSelectorPane.
	SEE END OF DOCUMENT FOR DETAILS. */

	var parameterWidget = this.parameterWidget;

	// SET SELECT CALLBACK
	var selectCallback = dojo.hitch(this, function(file, location, type, parameterWidget)
	{
        console.dir(parameterWidget);

		var newValue;
		if ( file != null && location != null )	newValue = location + "/" + file;
		else if ( location != null )	newValue = location;
		else if ( file != null )	newValue = file;

		parameterWidget.changeValue(parameterWidget.valueNode, parameterWidget.valueNode.innerHTML, newValue, type);
	});

	// SET SELECT CALLBACK
	var addCallback = dojo.hitch(this, function(file, location, type)
	{

		var newValue;
		if ( file != null && location != null )	newValue = location + "/" + file;
		else if ( location != null )	newValue = location;
		else if ( file != null )	newValue = file;

		parameterWidget.addValue(parameterWidget.valueNode, parameterWidget.valueNode.innerHTML, newValue, type);
	});


	this.fileManager = new plugins.workflow.FileManager(
	{
		paneId: "projects" + this.paneId,
		tabsNodeId: "tabs",
		selectCallback: selectCallback,
		addCallback: addCallback,
		parentWidget: this,
		core: this.core
	});

	this.core.fileManager = this.fileManager;
},

setRunStatus : function () {
// SET STATUS TAB IN INFO PANE BY INSTANTIATING RunStatus OBJECT

	if ( this.runStatus == null )
	{
		this.runStatus = new plugins.workflow.RunStatus(
		{
			attachNode : this.rightPane,
			parentWidget: this,
			core: this.core
		});

		this.core.runStatus = this.runStatus;
	}
},


setParameters : function () {
// SET DATA TAB IN INFO PANE BY INSTANTIATING Parameters OBJECT

	if ( this.parameterPane == null )
	{
		this.parameterPane = new plugins.workflow.Parameters(
		{
			attachNode : this.rightPane,
			parentWidget: this,
			core: this.core

		});

		this.core.parameters = this.parameterPane;
	}
},


setApps : function() {
// SET APPLICATIONS DND DRAG SOURCE  - A SERIES OF TITLE PANES CONTAINING APPLICATIONS
// GROUPED BY APPLICATION TYPE - AND CONTEXT MENU FOR APPLICATIONS

	// CREATE AN Apps OBJECT AND SLOT IT INTO this.apps
	this.apps = new plugins.workflow.Apps(
		{
			tabContainer : this.leftPane,
			parentWidget : this,
			core: this.core
		}
	);

	this.core.apps = this.apps;

},	// 	setApps


setSharedApps : function() {
// SET APPLICATIONS DND DRAG SOURCE  - A SERIES OF TITLE PANES CONTAINING APPLICATIONS
// GROUPED BY APPLICATION TYPE - AND CONTEXT MENU FOR APPLICATIONS

	// CREATE AN Apps OBJECT AND SLOT IT INTO this.apps
	this.sharedApps = new plugins.workflow.SharedApps(
		{
			tabContainer : this.leftPane,
			parentWidget : this,
			core: this.core
		}
	);

	this.core.sharedApps = this.sharedApps;
},	// 	setSharedApps

setStages : function() {
// SET WORKFLOWS DND DRAG SOURCE  - A SERIES OF TITLE PANES CONTAINING APPLICATIONS,
// GROUPED BY APPLICATION TYPE - AND ITS CONTEXT MENU

	// CREATE AN Apps OBJECT AND SLOT IT INTO this.stages
	this.stages = new plugins.workflow.Stages(
		{
			tabContainer : this.middlePane,
			parentWidget : this,
			core: this.core
		}
	);

	this.core.stages = this.stages;
},	// 	setStages

setHistory : function () {
// SET INFO PANE BY INSTANTIATING HISTORY OBJECT

	if ( this.historyPane == null )
	{
		this.historyPane = new plugins.workflow.History(
		{
			attachNode : this.middlePane,
			parentWidget: this,
			core: this.core
		});
	}

	this.core.historyPane = this.historyPane;
},

setSharedProjects : function () {
// SET SHARED WORKFLOWS TAB

	if ( this.sharedPane == null )
	{
		this.sharedPane = new plugins.workflow.SharedProjects(
		{
			tabContainer: this.middlePane,
			parentWidget: this,
			core: this.core
		});
		this.core.sharedProjects = this.sharedProjects;
	}
}


}); // end of plugins.workflow.Workflow

