dojo.provide("plugins.dnd.Target");

dojo.require("dojo.dnd.Source");

dojo.declare("plugins.dnd.Target",
	dojo.dnd.Target,
{

	/////}

// summary: a Target object, which can be used as a DnD target

// contextMenu: dijit.Menu OBJECT. CONTEXT MENU FOR ALL NODES
contextMenu : null,

// parentWidget: plugins.workflow.Workflows OBJECT.	
parentWidget : null,

constructor: function(node, params){


	this.contextMenu = params.contextMenu;

	this.parentWidget = params.parentWidget;

	// summary:
	//		a constructor of the Target --- see the `dojo.dnd.Source.constructor` for details
	this.isSource = false;
	dojo.removeClass(this.node, "dojoDndSource");
},

// RESET this.droppingApp TO FALSE
resetDropping : function () {
	this.droppingApp = false;
},

onDndDrop : function(source, nodes, copy) {
// OVERRIDE onDndDrop: USE dojo.connect TO ADD EVENT TO NEW ITEM
// DROPPED NODES MUST HAVE AN application SLOT

	//console.dir(this);
	////console.dir(source);

	// SANITY
	if ( nodes[0].application == null )	return;

	// RETURN IF DROP FLAG IS SET
	if ( this.droppingApp == true )
	{
		return;
	}

	// SET DROP FLAG
	this.droppingApp = true;

	// summary: topic event processor for /dnd/drop,
	// called to finish the DnD operation break box
	var newNode;
	do
	{ 
		if ( this.containerState != "Over" )
		{
			break;
		}

		var oldCreator = this._normalizedCreator;

		// transferring nodes from the source to the target
		if ( this != source )
		{
			// CLONE THE DROPPED NODE AND ADD THE
			// CLONE TO THE DROP TARGET
			this._normalizedCreator = function(node, hint)
			{
				var t = source.getItem(node.id);
				var n = node.cloneNode(true);
				n.parentWidget = node.parentWidget;
				n.id = dojo.dnd.getUniqueId();
				return {node: n, data: t.data, type: t.type};
			};
		}  

		// CLEAN UP - REMOVE SELECTION AND ANCHOR STYLE
		this._removeSelection();
		if ( this != source )
		{
			this._removeAnchor();
		}

		if ( this != source && !copy && !this.creator )
		{
			source.selectNone();
		}

		// INSERT DROPPED NODE INTO DROP TARGET
		this.insertNodes(true, nodes, this.before, this.current);
		this.sync();

		// COMPLETE THE NODE COPY:
		//
		// 1. TRANSFER THE METADATA FROM THE DROPPED NODE TO
		// THE CLONED NODE.
		// 
		// 2. INCREMENT BY ONE THE number OF THE NODES AFTER
		// THE INSERTION POINT OF THE NEW NODE.
		var belowInsertedNode = false;
		var allNodes = this.getAllNodes();

		//var allNodes = this.parentWidget.dropTargetContainer.childNodes;

		var thisObject = this;
		dojo.forEach(allNodes, function(node, i)
		{
			if ( node.application == null )
			{
				// CLONE THIS OTHERWISE GET AN INTERESTING ERROR
				// WHEN DUPLICATE COPIES OF THE SAME APPLICATION
				// ARE DROPPED (SHARING THE SAME application OBJECT)
				node.application = dojo.clone(nodes[0].application);

				// ADD appname TO APPLICATION
				node.application.appname = node.application.name;

				// ADD NUMBER TO APPLICATION
				// CAST number TO STRING FOR LATER SORTING
				node.application.appnumber = (i + 1).toString();
				node.application.number = (i + 1).toString();
				node.number = (i + 1).toString();

				// SET DEFAULT CLUSTER IS EMPTY
				if ( node.application["cluster"] == null )	node.application["cluster"] = '';
				// ADD PROJECT AND WORKFLOW TO node's APPLICATION
				node.application.project = thisObject.parentWidget.getProject();
				node.application.workflow = thisObject.parentWidget.getWorkflow();

				// SET USERNAME
				node.application.username = Agua.cookie('username');

				// INSTANTIATE SOURCE ROW 
				var stageRow = new plugins.workflow.StageRow(node.application);

				// SET CORE WORKFLOW OBJECTS
				stageRow.core = thisObject.core;

				stageRow.workflowWidget = thisObject.parentWidget.parentWidget;

				// CLEAR NODE CONTENT
				node.innerHTML = '';

				// APPEND stageRow WIDGET TO NODE
				node.appendChild(stageRow.domNode);

				// ADD CONTEXT MENU TO NODE
				thisObject.contextMenu.bind(node);

				// SET stageRow AS node.parentWidget ATTRIBUTE FOR ACCESS LATER:
				// --- (ALSO ADDED this.name.parentWidget = this IN StageRow.startup())
				//
				// 1. WHEN CALLING Workflow.loadParametersPane SO THAT THE CORRECT
				// StageRow HAS ITS validInputs SET ACCORDING TO THE OUTCOME
				// OF Workflow.loadParametersPane
				//
				// 2. FOR RESETTING OF number ON REMOVAL OR INSERTION OF NODES
				//
				// REM: remove ONCLICK BUBBLES ON stageRow.name NODE RATHER THAN ON node. 
				// I.E., CONTRARY TO DESIRED, thisObject.name IS THE TARGET INSTEAD OF THE node.
				node.parentWidget = stageRow;

				//NB: NOT THIS: node.parentWidget = dojo.clone(nodes[0].parentWidget);

				// INSERT STAGE INTO thisObject.stage AND ITS STAGE PARAMETERS
				// INTO Agua.stages AND Agua.stageparameters
				var insertOk = Agua.insertStage(node.application);
				if ( ! insertOk )
				{
					// UNSET droppingApp FLAG
					//thisObject.droppingApp = false;
					setTimeout(thisObject.resetDropping, 1000);
					return;
				}

				// ADD ONCLICK TO LOAD APPLICATION INFO
				node.onclick = function(e)
				{
					thisObject.parentWidget.loadParametersPane(node, null);
				}

				// ADD CONTEXT MENU
				thisObject.contextMenu.bind(node);

				// SET THE DEFAULT CHAINED VALUES FOR INPUTS AND OUTPUTS FOR THE
				// APPLICATION BASED ON THOSE OF THE PREVIOUS APPLICATIONS
				thisObject.parentWidget.getChainedValues(node);


//*********** DEBUG ****************//
//*********** DEBUG ****************//

				// CHECK IF ALL STAGES ARE VALID.
				// NB: TIMEOUT DELAY REQUIRED TO ALLOW TIME
				// FOR COMPLETION OF NODE REMOVAL/ADDITION 
				// AND OF checkFiles IN INFOPANE
				//setTimeout(function(thisObj) { thisObj.checkValidStages(); }, 500, thisObject.parentWidget);


//*********** DEBUG ****************//
//*********** DEBUG ****************//


				// SET belowInsertedNode FLAG TO TRUE
				belowInsertedNode = true;

				// UNSET droppingApp FLAG
				thisObject.droppingApp = false;

				// SET NEW NODE FOR LOAD INFO PANE LATER
				newNode = node;


				return;
			}

			// IF WE ARE BELOW THE INSERTED NODE, CHAIN THE STAGE
			if ( belowInsertedNode == true )
			{

				var force = true;
				thisObject.parentWidget.workflowIO.chainStage(node.application, force);
			}
		});

		// UNSET droppingApp FLAG
		thisObject.droppingApp = false;

		thisObject._normalizedCreator = oldCreator;
	}
	while(false);
	// end of 'do'

	// SET THE APPLICATION.NUMBER AND .APPNUMBER FOR EACH NODE AND ITS WIDGET
	this.parentWidget.resetNumbers();


// ******************* DISABLED FOR DEBUGGING ***********************
// ******************* DISABLED FOR DEBUGGING ***********************


	// SET INFO PANE FOR DROPPED NODE
	this.parentWidget.loadParametersPane(newNode);


// ******************* DISABLED FOR DEBUGGING ***********************
// ******************* DISABLED FOR DEBUGGING ***********************


	this.onDndCancel();

}	// OVERRIDE onDndDrop TO USE dojo.connect TO ADD EVENT TO NEW ITEM







});
