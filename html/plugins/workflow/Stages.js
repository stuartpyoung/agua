dojo.provide("plugins.workflow.Stages");

/* 
	INPUT VALIDITY CHECKING PROCESS ON LOAD NEW WORKFLOW:

	1. setDropTarget: load stages

	2. setDropTarget: CALL -> first stage)

		   load multiple ParameterRows, each checks isValid (async xhr request for each file)
			   CALL -> Agua.setParameterValidity(boolean) to set stageParameter.isValid

	3. setDropTarget: (concurrently with 2.) CALL -> updateValidity()

		   check isValid for each StageRow:

				CALL-> checkValidParameters (async batch xhr request for multiple files)

				   CALL -> Agua.getParameterValidity(), and if empty check input and then

						CALL -> Agua.setParameterValidity(boolean) to set stageParameter.isValid

*/

//dojo.require("dijit.dijit"); // optimize: load dijit layer
dojo.require("dijit.form.Button");
dojo.require("dijit.form.TextBox");
dojo.require("dijit.form.Textarea");
dojo.require("dojo.parser");
dojo.require("dojo.dnd.Source");

/* NOTES EDITOR
dojo.require("dijit.Editor");
dojo.require("dijit.form.DateTextBox");
dojo.require("dijit.form.Textarea");

dojo.require("dijit.form.TextBox");
dojo.require("dijit.form.ValidationTextBox");
dojo.require("dijit.form.NumberTextBox");
dojo.require("dijit.form.CurrencyTextBox");
dojo.require("dojo.currency");
dojo.require("dijit.Dialog");
*/

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

// TIMER
dojo.require("dojox.timing");

// TOOLTIP
dojo.require("dijit.Tooltip");

// TOOLTIP DIALOGUE
dojo.require("dijit.Dialog");
dojo.require("dijit.form.Textarea");
dojo.require("dijit.form.Button");

// WIDGETS IN TEMPLATE
dojo.require("dijit.layout.SplitContainer");
dojo.require("dijit.layout.ContentPane");

// INHERITS
dojo.require("plugins.core.Common");

// HAS A
dojo.require("plugins.workflow.StageRow");
dojo.require("plugins.workflow.StagesMenu");
dojo.require("plugins.workflow.IO");
dojo.require("plugins.core.ComboBox");
dojo.require("plugins.core.Confirm");
dojo.require("plugins.dnd.Target");

// INPUT DIALOG
dojo.require("plugins.core.InteractiveDialog");
dojo.require("plugins.core.SelectiveDialog");

dojo.declare(
    "plugins.workflow.Stages",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "workflow/templates/stages.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

//addingApp STATE
addingApp : false,

// OR USE @import IN HTML TEMPLATE
cssFiles : [ dojo.moduleUrl("plugins") + "/workflow/css/stages.css"],

// PARENT WIDGET
parentWidget : null,

// TAB CONTAINER
tabContainer : null,

// CONTEXT MENU
contextMenu : null,

// CORE WORKFLOW OBJECTS
core : null,

// PREVENT DOUBLE CALL ON LOAD
workflowLoaded : null,
dropTargetLoaded : null,

/////}

constructor : function(args) {

	// GET INFO FROM ARGS
	this.core = args.core;
	this.core.stages = this;

	this.parentWidget = args.parentWidget;
	this.tabContainer = args.tabContainer;

	// LOAD SORIA AND FILEPICKER CSS
	this.loadCSS();		
},

postCreate : function() {

	this.startup();
},


startup : function () {

	// COMPLETE CONSTRUCTION OF OBJECT
	this.inherited(arguments);	 

//return;

	// ADD TO TAB CONTAINER		
	this.tabContainer.addChild(this.mainTab);
	this.tabContainer.selectChild(this.mainTab);

	// SET INPUT/OUTPUT CHAINER
	this.setWorkflowIO();

	// SET SELECTIVE DIALOG FOR copyWorkflow	
	this.setSelectiveDialog();

	// SET INTERACTIVE DIALOG FOR copyProject
	this.setInteractiveDialog();

	// CREATE SOURCE MENU
	this.setContextMenu();

	// START CASCADE OF COMBO LOADING:
	// PROJECT COMBO, WORKFLOW COMBO, DROP TARGET
	this.setProjectCombo();

	// SET ONCLICK FOR PROJECT AND WORKFLOW BUTTONS
	this.setComboButtons();
},

setContextMenu : function () {
// GENERATE CONTEXT MENU


	this.contextMenu = new plugins.workflow.StagesMenu(
		{
			parentWidget: this,
			core: this.core
		}
	);
},


// DROP TARGET METHODS
setDroppingApp : function (value) {
// SET this.droppingApp


	this.droppingApp = value; 
},


setDropTarget : function (project, workflow) {
// SET DND DROP TARGET AND ITS CONTEXT MENU.
// ADD application OBJECT TO EACH NODE, CONTAINING
// THE STAGE INFORMATION INCLUDING number

	if ( this.dropTargetLoaded == null )
	{
		this.dropTargetLoaded = true;
		return;
	}

	// GET STAGES FOR THIS WORKFLOW
	var stages = Agua.getStagesByWorkflow(project, workflow);
	if ( stages == null )	stages = [];

	this._setDropTarget(stages);
},

updateDropTarget : function (workflow) {
// FIRED FROM onChange LISTENER IN workflowCombo.
// UPDATE THE STAGES IN THE DROP TARGET
	var project = this.getProject();

	this.setDropTarget(project, workflow);
},

_setDropTarget : function (stages) {

	//return;

	// EMPTY dropTargetContainer
	while ( this.dropTargetContainer.firstChild )
	{
		this.dropTargetContainer.removeChild(this.dropTargetContainer.firstChild);
	}

	// SORT STAGES
	stages = this.sortNumericHasharray(stages, "number");

	// GENERATE APPLICATIONS ARRAY FOR DRAG AND DROP
	// FROM WORKFLOW APPLICATIONS LIST
	var dataArray = new Array;
	for ( var i = 0; i < stages.length; i++ )
	{
		var hash = new Object;
		hash.data = stages[i].name;
		hash.type = [ "draggableItem" ];
		dataArray.push(hash);
	}

	// CREATE DROP TARGET
	if ( this.dropTarget == null )
	{
		this.dropTarget = new plugins.dnd.Target( this.dropTargetContainer,
			{
				accept: ["draggableItem"],
				contextMenu : this.contextMenu,
				parentWidget : this,
				core: this.core
			}
		);
	}

	// DELETE ALL NODES	(NEEDED TO CLEAN UP REMAINING NODES FROM PREVIOUS WORKFLOW)	
	//var selectedNodes = this.dropTarget.selectAll();
	//this.dropTarget.deleteSelectedNodes(selectedNodes);

	// INSERT DATA INTO DROP TARGET
	this.dropTarget.insertNodes(false, dataArray);

	// SET NODE CHARACTERISTICS - ONCLICK, CLASS, ETC.
	allNodes = this.dropTarget.getAllNodes();
	var thisObject = this;

	dojo.forEach(allNodes, function (node, i)
	{

		var nodeClass = dataArray[i].type;
		var applicationName = node.innerHTML;

		// ADD infoId TO NODE
		node.setAttribute('infoId', thisObject.infoId);

		// ADD application TO NODE
		node.application = stages[i];

		// DON'T SHOW DESCRIPTION ONCLICK FOR NOW
		// BECAUSE WE ALREADY HAVE AN ONCLICK
		// LISTENER FOR REFRESHING THE INFOPANE
		node.application.description = '';

		// NO DESCRIPTION OR NOTES FOR NOW IN STAGE
		stages[i].description = '';
		stages[i].notes = '';

		// NO INDIVIDUAL CLUSTER FOR STAGE AS YET
		if ( stages[i]["cluster"] == null )	stages[i]["cluster"] = '';

		// SET PARENT WIDGET
		//stages[i].parentWidget = this;

		// INSTANTIATE SOURCE ROW 
		var stageRow = new plugins.workflow.StageRow(stages[i]);
		stageRow.core = thisObject.core;

		// SET stageRow.parentWidget = Stages
		stageRow.parentWidget = thisObject;

		// CLEAR NODE CONTENT
		node.innerHTML = '';

		// APPEND TO NODE
		node.appendChild(stageRow.domNode);

		// SET stageRow AS node.parentWidget FOR LATER RESETTING OF
		// number ON REMOVAL OR INSERTION OF NODES
		//
		// REM: remove ONCLICK BUBBLES ON stageRow.name NODE RATHER THAN ON node. 
		// I.E., CONTRARY TO DESIRED, this.name IS THE TARGET INSTEAD OF THE node.
		//
		// ALSO ADDED this.name.parentWidget = this IN StageRow.startup()
		node.parentWidget = stageRow;

		// ADD CONTEXT MENU TO NODE
		thisObject.contextMenu.bind(node);

		// SHOW APPLICATION INFO WHEN CLICKED
		dojo.connect(node, "onclick",  dojo.hitch(thisObject, function(event)
			{

				event.stopPropagation();
				this.loadParametersPane(node);
			}
		));

	});	// END OF allNodes


	// 1. SET THE INFOPANE FOR THE FIRST STAGE:
	//		- CHECK VALIDITY OF PARAMETERS
	//		- CHANGE PARAMETER NODE STYLES ACCORDINGLY
	// 2. CHECK VALIDITY OF PARAMETERS FOR REMAINING STAGES
	if ( allNodes.length != 0 )
	{



//////////     DEBUG    /////////////		
//////////     DEBUG    /////////////		
//////////     DEBUG    /////////////		


		// Workflow.loadParametersPane
		//     --> CALLS Parameters.load
		//         --> CALLS StageRow.checkValidParameters
		// CARRIED OUT SYNCHRONOUSLY (I.E., WAITS TIL DONE)
		this.loadParametersPane(allNodes[0]);

		// CALL StageRow.checkValidParameters DIRECTLY FOR 
		// THE REMAINING STAGES
		// (I.E., IGNORING THE PARAMETER PANE)
		if ( allNodes.length >= 1 )
		{
			for ( var i = 1; i < allNodes.length; i++ )
			{
				var stageRow = allNodes[i].parentWidget;

				// DON'T FORCE IN CASE STAGE PARAMETER INFORMATION 
				// HAS ALREADY BEEN GENERATED EARLIER IN THIS SESSION
				var force = false;
				stageRow.checkValidParameters(force);
			}	
		}
	}

	// UPDATE THE VALIDITY OF Stages BASED ON VALIDITY OF StageRows	
	this.updateValidity();

	// CHECK ALL STAGES ARE VALID, WITH TIMEOUT FOR checkFiles IN INFOPANE	
//	setTimeout(function(thisObj) { thisObj.updateValidity(); }, 100, this);


	// CHECK RUN STATUS, JUST IN CASE THE WORKFLOW IS ALREADY RUNNING
	//this.runStatus();

}, // end of Stages._setDropTarget


resetNumbers : function () {

	var childNodes = this.dropTarget.getAllNodes();

	// RESETTING number IN ALL CHILDNODES
	for ( var i = 0; i < childNodes.length; i++ )
	{
		var node = childNodes[i];

		// GET WIDGET
		var widget = dijit.byNode(node.firstChild);
		//var widget = node.parentWidget;
		if ( widget == null )
		{
			widget = dijit.getEnclosingWidget(childNodes[i]);
		}
		node.application.number = (i + 1).toString();
		node.application.appnumber = (i + 1).toString();

		widget.setNumber(node.application.number);

	}
},


// CLUSTER METHODS
getCluster : function () {
	return this.clusterCombo.get('value');
},

toggleClusterCombo : function () {
	var disabled = this.clusterCombo.get('disabled');
	if ( disabled == true )
	{
		this.clusterCombo.set('disabled', false);
		dojo.addClass(this.clusterCombo.comboNode, 'enabled')
	}
	else
	{
		dojo.removeClass(this.clusterCombo.comboNode, 'enabled')
		this.clusterCombo.set('disabled', true);
	}
},

setClusterCombo : function () {
// POPULATE THE WORKFLOW COMBO BASED ON SELECT VALUE IN PROJECT COMBO

	var projectName = this.getProject();
	var workflowName = this.getWorkflow();

	var clusters = Agua.getClusters();

    var cluster = Agua.getClusterByWorkflow(projectName, workflowName);
	if ( cluster == null )	cluster = '';

	var clusterNames = this.hasharrayKeyToArray(clusters, "cluster");
	clusterNames.splice(0,0, '');	

	// DO data FOR store
	var data = {identifier: "name", items: []};
	for ( var i in clusterNames )
	{
		data.items[i] = { name: clusterNames[i]	};
	}

	// CREATE store
	var store = new dojo.data.ItemFileReadStore( { data: data } );

	// GET WORKFLOW COMBO
	if ( this.clusterCombo == null )
	{
		return;
	}

	this.clusterCombo.store = store;
	this.clusterCombo.startup();
	this.clusterCombo.set('value', cluster);			

	// SET CSS
	this.clusterCombo.popupClass = "workflow stages dijitReset dijitMenu";
	this.clusterCombo.wrapperClass = "stages dijitPopup";
	this.clusterCombo.itemHeight = 30;
},

checkEnterNodes : function (event) {

	if (event.keyCode == dojo.keys.ENTER)
	{
		document.body.focus();

		this.checkNodeNumbers();

		this.updateCluster();
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

updateCluster : function () {
//	SAVE A PARAMETER TO Agua.parameters AND TO REMOTE DATABASE


	if ( this.savingCluster == true )	return;
	this.savingCluster = true;

	var clusterObject = new Object;
	clusterObject.username = Agua.cookie('username');	
	clusterObject.sessionId = Agua.cookie('sessionId');	
	clusterObject.minnodes = this.minNodes.get('value');
	clusterObject.maxnodes = this.maxNodes.get('value');
	clusterObject.cluster = this.getCluster();	

	Agua.removeCluster(clusterObject); // not necessary as addCluster does delete first
	Agua.addCluster(clusterObject);

	this.savingCluster = false;

	// SAVE ON REMOTE DATABASE
	var url = Agua.cgiUrl + "/agua?";
	clusterObject.mode = "updateCluster";

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

}, // Clusters.updateCluster


saveClusterWorkflow : function () {
//	SAVE A PARAMETER TO Agua.parameters AND TO REMOTE DATABASE

	// AVOID DOUBLE CALL ON INSTANTIATION
	if ( this.clusterLoaded == null )
	{
		this.clusterLoaded = true;
		return;
	}

	if ( this.savingCluster == true )	return;
	this.savingCluster = true;

	var clusterObject = new Object;
	clusterObject.username = Agua.cookie('username');	
	clusterObject.sessionId = Agua.cookie('sessionId');	
	clusterObject.project = this.getProject();
	clusterObject.workflow = this.getWorkflow();
	clusterObject.cluster = this.getCluster();	

	Agua.removeCluster(clusterObject); // not necessary as addCluster does delete first
	Agua.addCluster(clusterObject);

	this.savingCluster = false;

	// SAVE ON REMOTE DATABASE
	var url = Agua.cgiUrl + "/agua?";
	clusterObject.mode = "saveClusterWorkflow";

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

}, // Clusters.saveClusterWorkflow


/////////////////		INFOPANE METHODS         
clearParameters : function () {
// CLEAR INFO PANE

	if ( this.parentWidget.parameterPane == null )
	{
		return;
	}

	this.parentWidget.parameterPane.clear();
},


loadParametersPane : function (node) {
// LOAD DATA INTO INFO PANE FROM THE APPLICATION ASSOCIATED WITH THIS NODE
// OVERLOAD THIS TO PASS ADDITIONAL ARGUMENTS TO Parameters.load()


	// WARN AND QUIT IF NO NODE PASSED, E.G., IF WORKFLOW HAS NO STAGES
	if ( node == null )
	{
		return;
	}

	if ( this.parentWidget.parameterPane != null )
	{
		this.parentWidget.parameterPane.load(node);
	}
},


/////////////////		WorkflowIO METHODS      
setWorkflowIO : function () {
// INITIATE this.workflowIO OBJECT

	if ( this.workflowIO == null )
	{
		this.workflowIO = new plugins.workflow.IO(
			{
				parentWidget: this,
				core: this.core
			}
		);
	}
},


getChainedValues : function (node) {
// SET THE input, resource AND output PARAMETERS OF THIS STAGE USING
// ANY CORRESPONDING PARAMETERS IN THE PRECEDING STAGE


	// CHECK node.application AND node.application.number ARE DEFINED
	// GET THE INDEX OF THIS APPLICATION
	if ( node.application == null )
	{
		return;
	}

	// GET THE INDEX OF THIS APPLICATION
	if ( node.application.number == null )
	{
		this.quit("Stages.getChainedValues    node.application.number is null. Returning.");
		return;
	}

	// CHANGE THE STAGE PARAMETERS FOR THIS APPLICATION
	// IF THE args FIELD IS NOT NULL (ALSO params AND paramFunction)
	var force = true;
	this.workflowIO.chainStage(node.application, force);
},

chainStages : function (force) {
// CHAIN THE INPUTS/OUTPUTS OF ALL APPLICATIONS IN THE WORKFLOW
// NB: NOT USED YET, JUST IN CASE


	var nodes = this.dropTarget.getAllNodes();

	for ( var i = 0; i < nodes.length; i++ )
	{
		this.workflowIO.chainInputs(nodes[i].application, force);
	} 
},


/////////////////		PROJECT METHODS      
getProject : function () {

	return this.projectCombo.get('value');
},

setProjectCombo : function (project, workflow) {
// POPULATE THE PROJECT COMBO AND THEN RELOAD THE WORKFLOW COMBO

	//
	this.inherited(arguments);

	// SET CSS
	this.projectCombo.popupClass = "workflow stages dijitReset dijitMenu";
	this.projectCombo.wrapperClass = "stages dijitPopup";
	this.projectCombo.itemHeight = 30;

	// SET NEW PROJECT LISTENER
	var thisObject = this;
	this.projectCombo._onKeyPress = function(evt){

		// summary: handles keyboard events
		var key = evt.charOrCode;

		if ( key == 13 )
		{
			thisObject.projectCombo._hideResultList();

			var sourceProject = thisObject.projectCombo.getValue();

			if ( Agua.isProject(sourceProject) == false )
			{
				// CLEAR THE INFO PANE
				thisObject.parentWidget.clearParameters();

				// ADD THE PROJECT
				Agua.addProject(sourceProject);

				// RESET THE PROJECT COMBO
				thisObject.setProjectCombo(sourceProject);
			}

			if ( thisObject.projectCombo._popupWidget != null )
			{
				thisObject.projectCombo._showResultList();
			}
		}
	};


	// DOJO.CONNECT TO CHANGE THE workflowCombo
	var thisObject = this;
	dojo.connect(this.projectCombo, "onChange", dojo.hitch(this, function(event) {

		// CLEAR THE INFO PANE
		thisObject.clearParameters();

		// SET WORKFLOW COMBO
		thisObject.project = event;
		thisObject.setWorkflowCombo(event);
	}));


	if ( project == null )
	{
		project = this.projectCombo.getValue();
	}

	// RESET THE WORKFLOW COMBO
	this.setWorkflowCombo(project, workflow);
},

deleteProject : function (event) {
// DELETE A PROJECT AFTER ONCLICK 'DELETE' BUTTON BESIDE project COMBOBOX


	// SET this.doingDelete OR EXIT IF BUSY
	if ( this.doingDelete == true )
	{
		return;
	}
	this.doingDelete = true;

	if ( ! Agua.stages )
	{
		return;
	}

	// GET DELETED PROJECT NAME OR QUIT IF EMPTY
	var sourceProject = this.projectCombo.getValue();
	if ( sourceProject == null || ! sourceProject )
	{
		this.doingDelete = false;
		return;
	}

	var projectObject = {
		name: sourceProject,
		workflow: null
	}
	// DELETE FROM AGUA
	Agua.removeProject(projectObject);

	// SET ARGS FOR CONFIRM DELETE
	var args = new Object;
	args.project = sourceProject;
	args.workflow = null;

	// DO CONFIRM DELETE
	this.confirmDelete(args);

	// UNSET this.doingDelete
	this.doingDelete = false;
},


/////////////////		WORKFLOW METHODS      
getWorkflow : function () {
	return this.workflowCombo.get('value');
},

getWorkflowNumber : function () {
	var project = this.getProject();
	var workflow = this.getWorkflow();
	if ( project == null || workflow == null )
	{
		return;
	}

	return Agua.getWorkflowNumber(project, workflow);
},

setWorkflowCombo : function (projectName, workflowName) {
// POPULATE THE WORKFLOW COMBO BASED ON SELECT VALUE IN PROJECT COMBO


	// AVOID DOUBLE CALL ON INSTANTIATION
	if ( this.workflowLoaded == null )
	{
		this.workflowLoaded = true;
		return;
	}

	// POPULATE THE WORKFLOW COMBO AND SET FIRST VALUE TO
	// workflowName OR THE FIRST WORKFLOW IF workflowName NOT DEFINED
	this.inherited(arguments);


	// SET CSS
	this.workflowCombo.popupClass = "workflow stages dijitReset dijitMenu";
	this.workflowCombo.wrapperClass = "stages dijitPopup";
	this.workflowCombo.itemHeight = 30;

	// SET NEW PROJECT LISTENER
	var thisObject = this;
	this.workflowCombo._onKeyPress = function(evt){

		// summary: handles keyboard events
		var key = evt.charOrCode;			
		if ( key == 13 )
		{
			thisObject.workflowCombo._hideResultList();

			var projectName = thisObject.getProject();
			var workflowName = thisObject.getWorkflow()

			// STOP PROPAGATION
			evt.stopPropagation();

			var isWorkflow = Agua.isWorkflow(projectName, workflowName);
			if ( isWorkflow == false )
			{
				Agua.addWorkflow(projectName, workflowName);
				thisObject.setWorkflowCombo(projectName, workflowName);
			}

			if ( thisObject.workflowCombo._popupWidget != null )
			{
				thisObject.workflowCombo._showResultList();
			}
		}
	};

	// DOJO.CONNECT TO POPULATE APPLICATIONS IN DROP TARGET
	// WHICH THEN POPULATES THE INFO PANE 
	var thisObject = this;
	dojo.connect(thisObject.workflowCombo, "onChange", dojo.hitch(this, function(event) {
		var workflowName = event;
		var projectName = thisObject.getProject();

		thisObject.setDropTarget(projectName, workflowName);
	}));

	// SET DROP TARGET (LOAD MIDDLE PANE, BOTTOM)
	if ( workflowName == null )
	{
		workflowName = this.workflowCombo.getValue();
	}

	this.setClusterCombo(projectName, workflowName);

	this.setDropTarget(projectName, workflowName);
},


newWorkflow : function (sourceProject, workflowName) {
// CREATE A NEW WORKFLOW ON TRIGGER this.workflowCombo._onKeyPress ENTER


	if ( this.doingNewWorkflow == true )
	{
		return;
	}

	// SET this.doingNewWorkflow
	this.doingNewWorkflow = true;

	// SEND TO SERVER
	Agua.addWorkflow(sourceProject, workflowName);

	// UNSET this.doingNewWorkflow
	this.doingNewWorkflow = false;

	// RESET THE WORKFLOW COMBO
	this.setWorkflowCombo(sourceProject, workflowName);

	// SEND TO SERVER
	Agua.addProjectWorkflow(sourceProject, workflowName);
},




setComboButtons : function () {
// SET ONLICK LISTENERS FOR PROJECT AND WORKFLOW DELETE BUTTONS

	// SET download BUTTON ONCLICK TO OPEN FILE MANAGER
	var thisObject = this;

	dojo.connect(this.deleteProjectButton, "onclick", function(event)
	{
		thisObject.deleteProject(event);
	});

	dojo.connect(this.deleteWorkflowButton, "onclick", function(event)
	{
		thisObject.deleteWorkflow(event);
	});
},



deleteWorkflow : function (event) {
// DELETE A PROJECT AFTER ONCLICK 'DELETE' BUTTON BESIDE project COMBOBOX


	// SET this.doingDelete OR EXIT IF BUSY
	if ( this.doingDelete == true )
	{
		return;
	}
	this.doingDelete = true;

	// GET DELETED PROJECT NAME OR QUIT IF EMPTY
	var project = this.projectCombo.getValue();
	if ( project == null || ! project )
	{
		this.doingDelete = false;
		return;
	}

	// GET DELETED WORKFLOW NAME OR QUIT IF EMPTY
	var workflow = this.workflowCombo.getValue();
	if ( workflow == null || ! workflow )
	{
		this.doingDelete = false;
		return;
	}

	// SET ARGS FOR CONFIRM DELETE
	var args = new Object;
	args.project = project;
	args.workflow = workflow;

	// DO CONFIRM DELETE
	this.confirmDelete(args);

	// UNSET this.doingDelete
	this.doingDelete = false;
},








/////////////////		STAGE METHODS         
updateStageNumber : function (stageObject, previousNumber) {
// UPDATE THE number OF A STAGE IN this.stages
// AND ON THE REMOTE SERVER


	// REMOVE FROM Agua DATA
	var addOk = Agua.updateStageNumber(stageObject, previousNumber);
	if ( ! addOk )
	{
		return;
	}
},


/////////////////		RUN STATUS METHODS         
updateValidity : function () {
// CHECK ALL STAGE INPUTS ARE VALID, ADJUST 'RUN' BUTTON CSS ACCORDINGLY


	var stageNodes = this.dropTarget.getAllNodes();

	this.isValid = true;

	for ( var i = 0; i < stageNodes.length; i++ )
	{
		var stageRow = stageNodes[i].parentWidget;
		//var stageRow = dijit.getEnclosingWidget(stageNodes[i]);
		if ( stageRow == null )
			stageRow = dijit.byNode(stageNodes[i].firstChild);

		if ( stageRow == null )
		{
			this.isValid = false;
			return;
		}


		if ( stageRow.isValid == false || stageRow.isValid == null )
			this.isValid = false;
	}	

	if ( this.isValid == true )	this.enableRunButton();
	else this.disableRunButton();
},


enableRunButton : function () {
// ENABLE RUN BUTTON - ADD ONCLICK AND REMOVE invalid CSS


	// GET RUN BUTTON AND TITLE NODE		
	var node = this.runButton;

	// ADD enabled CSS
	dojo.removeClass(node, 'runButtonDisabled');
	dojo.addClass(node, 'runButtonEnabled');

	// REMOVE 'RUN' ONCLICK
	if ( node.onclickListener != null )
	{
		dojo.disconnect(node.onclickListener);
	}

	// SET 'RUN' ONCLICK
	var thisObject = this;
	node.onclickListener = dojo.connect( node, "onclick", function(event)
	{
		var runner = new Object;
		runner.project = thisObject.getProject();
		runner.workflow = thisObject.getWorkflow();
		runner.cluster = thisObject.getCluster();
		runner.workflownumber = thisObject.getWorkflowNumber();
		runner.number = 1;
		runner.childNodes = thisObject.dropTarget.getAllNodes();

		// RUN ALL STAGES IN THE WORKFLOW (ASSUMES ALL STAGES ARE VALID)

		thisObject.parentWidget.runStatus.startRun(runner);
	});
},


disableRunButton : function () {
// DISABLE RUN BUTTON - REMOVE ONCLICK AND ADD invalid CSS


	// GET RUN BUTTON AND TITLE NODE		
	var node = this.runButton;

	// REMOVE enabled CSS
	dojo.removeClass(node, 'runButtonEnabled');
	dojo.addClass(node, 'runButtonDisabled');


	// REMOVE 'RUN' ONCLICK
	if ( node.onclickListener != null )
	{
		dojo.disconnect(node.onclickListener);
	}

},




/////////////////		CONFIRM DELETE
commitDelete : function (args) {
// DELETE THE WORKFLOW/PROJECT AND UPDATE THE WORKFLOW COMBO BOX

	if ( args.project == null )	return;

	// DELETE THE WORKFLOW AND UPDATE THE WORKFLOW COMBO BOX
	if ( args.workflow != null )
	{
		Agua.removeWorkflow({ project: args.project, name: args.workflow});

		// CLEAR THE INFO PANE
		this.clearParameters();

		//  RESET THE WORKFLOW COMBO
		var sourceProject = this.projectCombo.getValue();
		this.setWorkflowCombo(sourceProject);
	}

	// DELETE THE PROJECT AND UPDATE THE PROJECT COMBO BOX
	else
	{
		Agua.removeProject({ name: args.project });

		//  RESET THE WORKFLOW COMBO
		this.setProjectCombo();
	}			

	Agua.toastMessage("Deleted workflow '" + args.workflow + "' in project '" + args.project + "'");
},



confirmDelete : function (args) {

	// SET CALLBACKS
	var thisObject = this;
	var yesCallback = function()
	{
		thisObject.commitDelete(args);
	};
	var noCallback = function(){};

	// SET TITLE
	var title = "Delete project: " + args.project + "?";
	if ( args.workflow != null )
		title = "Delete workflow: " + args.workflow + "?";

	// SET MESSAGE
	var message = "All stages and data will be destroyed<br><span style='color: #222;'>Click 'Yes' to delete or 'No' to cancel</span>";
	if ( args.workflow != null )
		message = "All data will be destroyed<br><span style='color: #222;'>Click 'Yes' to delete or 'No' to cancel</span>";


	// THIS DOESN'T WORK YET: CALLBACKS CHANGE OK BUT NODE DISPLAYED VALUES DON'T
	// SO GO WITH REINSTANTIATION FOR NOW
	//
	//// IF NOT EXISTS, INSTANTIATE WIDGET CONTAINING CONFIRMATION DIALOGUE POPUP
	//if ( this.confirm != null )
	//{
	//	this.confirm.messageNode.innerHTML 	= message;
	//	this.confirm.yesCallback	= yesCallback;
	//	this.confirm.noCallback		= noCallback;
	//	this.confirm.dialog.titleNode.value	= title;
	//	this.confirm.dialog.reset();
	//	this.confirm.dialog.show();
	//}
	//
	//else
	//{
	//	// OTHERWISE, LOAD THE NEW VALUES AND SHOW THE DIALOGUE
	//	this.confirm = new plugins.core.Confirm({
	//		parentWidget : this,
	//		title: title,
	//		message : message,
	//		yesCallback : yesCallback,
	//		noCallback : noCallback
	//	});
	//}	



	// IF NOT EXISTS, INSTANTIATE WIDGET CONTAINING CONFIRMATION DIALOGUE POPUP
	if ( this.confirm != null ) 	this.confirm.destroy();

	// LOAD THE NEW VALUES AND SHOW THE DIALOGUE
	this.confirm = new plugins.core.Confirm({
		parentWidget : this,
		title: title,
		message : message,
		yesCallback : yesCallback,
		noCallback : noCallback
	});
	this.confirm.show();
},


setSelectiveDialog : function () {
	var enterCallback = function (){};
	var cancelCallback = function (){};
	var title = "";
	var message = "";

	this.selectiveDialog = new plugins.core.SelectiveDialog(
		{
			title 				:	title,
			message 			:	message,
			parentWidget 		:	this,
			enterCallback 		:	enterCallback,
			cancelCallback 		:	cancelCallback
		}			
	);
},

loadSelectiveDialog : function (title, message, comboValues, inputMessage, comboMessage, checkboxMessage, enterCallback, cancelCallback) {


	this.selectiveDialog.load(
		{
			title 				:	title,
			message 			:	message,
			comboValues 		:	comboValues,
			inputMessage 		:	inputMessage,
			comboMessage 		:	comboMessage,
			checkboxMessage		:	checkboxMessage,
			parentWidget 		:	this,
			enterCallback 		:	enterCallback,
			cancelCallback 		:	cancelCallback
		}			
	);
},



setInteractiveDialog : function () {
	var enterCallback = function (){};
	var cancelCallback = function (){};
	var title = "";
	var message = "";

	this.interactiveDialog = new plugins.core.InteractiveDialog(
		{
			title 				:	title,
			message 			:	message,
			parentWidget 		:	this,
			enterCallback 		:	enterCallback,
			cancelCallback 		:	cancelCallback
		}			
	);
},

loadInteractiveDialog : function (title, message, enterCallback, cancelCallback, checkboxMessage) {

	this.interactiveDialog.load(
		{
			title 				:	title,
			message 			:	message,
			checkboxMessage 	:	checkboxMessage,
			enterCallback 		:	enterCallback,
			cancelCallback 		:	cancelCallback
		}			
	);
},


/////////////////		COPY WORKFLOW / PROJECT
copyWorkflow : function () {
// DISPLAY A 'Copy Workflow' DIALOG THAT ALLOWS THE USER TO SELECT 
// THE DESTINATION PROJECT AND THE NAME OF THE NEW WORKFLOW


	// SET TITLE AND MESSAGE
	var sourceProject = this.projectCombo.get('value');
	var sourceWorkflow = this.workflowCombo.get('value');



	// SET CALLBACKS
	var cancelCallback = function (){
	};
	var thisObject = this;

	var enterCallback = dojo.hitch(this, function (targetProject, targetWorkflow, copyFiles, dialogWidget)
		{

			// SET BUTTON LABELS
			var enterLabel = "Copy";
			var cancelLabel = "Cancel";

			// SANITY CHECK
			if ( targetWorkflow == null || targetWorkflow == '' )	return;
			targetWorkflow = targetWorkflow.replace(/\s+/, '');

			// QUIT IF WORKFLOW EXISTS ALREADY
			if ( Agua.isWorkflow(targetProject, targetWorkflow) == true )
			{

				dialogWidget.messageNode.innerHTML = "/" + targetWorkflow + "' already exists in '" + targetProject + "'";
				return;
			}
			else {

				dialogWidget.messageNode.innerHTML = "Creating workflow";
				dialogWidget.close();
			}

			thisObject._copyWorkflow(sourceProject, sourceWorkflow, targetWorkflow, targetProject, copyFiles);
		}
	);		

	// SHOW THE DIALOG
	this.selectiveDialog.load(
		{
			title 				:	"Copy WorkflowXXX",
			message 			:	"Source: '" + sourceProject + ":" + sourceWorkflow + "'",
			comboValues 		:	Agua.getProjectNames(),
			inputMessage 		:	"Workflow",
			comboMessage 		:	"Project",
			checkboxMessage		:	"Copy files",
			parentWidget 		:	this,
			enterCallback 		:	enterCallback,
			cancelCallback 		:	cancelCallback,
			enterLabel			:	"Copy",
			cancelLabel			:	"Cancel"
		}			
	);

},

_copyWorkflow : function (sourceProject, sourceWorkflow, targetProject, targetWorkflow, copyFiles) {

	var username = Agua.cookie('username');
	// ADD PROJECT
	Agua.copyWorkflow(username, sourceProject, sourceWorkflow, username, targetProject, targetWorkflow, copyFiles);

	this.setWorkflowCombo(targetProject, targetWorkflow);
},

copyProject : function () {
// ADD A NEW WORKFLOW USING A DIALOG BOX FOR WORKFLOW NAME INPUT


	var sourceProject = this.projectCombo.get('value');

	// SET CALLBACKS
	var cancelCallback = function (){
	};

	var thisObject = this;
	var enterCallback = dojo.hitch(this, function (targetProject, copyFiles, interactiveDialog)
		{

			// SANITY CHECK
			if ( targetProject == null || targetProject == '' )	return;
			targetProject = targetProject.replace(/\s+/, '');

			// QUIT IF WORKFLOW EXISTS ALREADY
			if ( Agua.isProject(targetProject) == true )
			{

				interactiveDialog.messageNode.innerHTML = "Project name already exists";
				return;
			}
			else {

				interactiveDialog.messageNode.innerHTML = "Creating project";
				interactiveDialog.close();
			}

			thisObject._copyProject(sourceProject, targetProject, copyFiles);
		}
	);	

	// SHOW THE DIALOG
	this.interactiveDialog.load(
		{
			title 				:	"Copy Project",
			message 			:	"Please enter project name",
			checkboxMessage 	:	"Copy files",
			enterCallback 		:	enterCallback,
			cancelCallback 		:	cancelCallback
		}			
	);
},

_copyProject : function (sourceProject, targetProject, copyFiles) {

	var username = Agua.cookie('username');
	// ADD PROJECT
	Agua.copyProject(username, sourceProject, username, targetProject, copyFiles);

	this.setProjectCombo(targetProject);
}






}); // plugins.workflow.Stages
