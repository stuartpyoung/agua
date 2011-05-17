// OBJECT:  Workflow
// PURPOSE: CREATE A NEW TAB PANE 

dojo.declare( "Main", null,
{
	// PANE IDENTIFIER
	paneId : '',

	// ARRAY OF WORKFLOW APPLICATIONS
	workflowApplications : [],

	// CONSTRUCTOR	
	constructor : function(args) {

		if ( ! args.id )
		{
			return;
		}

		this.paneId = args.id;
	},

	// ADD PROGRAMMATIC CONTEXT MENU
	createMenu : function ()
	{
		var dynamicMenu = new dijit.Menu( { id: 'dynamicMenuPopup'} );

		// ADD MENU TITLE
		dynamicMenu.addChild(new dijit.MenuItem( { label:"Application Menu", disabled:false} ));
		dynamicMenu.addChild(new dijit.MenuSeparator());

		//// ONE OF FOUR WAYS TO DO MENU CALLBACK WITH ACCESS TO THE MENU ITEM AND THE CURRENT TARGET 	
		// 4. dojo.connect CALL
		//	REQUIRES:
		//		ADDED menu.currentTarget SLOT TO dijit.menu
		var mItem1 = new dijit.MenuItem(
			{
				id: "remove",
				label: "Remove",
				disabled:true 
			}
		);
		dynamicMenu.addChild(mItem1);
		dojo.connect(mItem1, "onClick", function()
			{
				var parentNode = dynamicMenu.currentTarget.parentNode;
				parentNode.removeChild(dynamicMenu.currentTarget);	
			}
		);

		// SEPARATOR
		dynamicMenu.addChild(new dijit.MenuSeparator());

		//	ADD run MENU ITEM
		var mItem2 = new dijit.MenuItem(
			{
				id: "run",
				label: "Run",
				disabled: true
			}
		);
		dynamicMenu.addChild(mItem2);	

		dojo.connect(mItem2, "onClick", function()
			{
				var currentTarget = dynamicMenu.currentTarget; 
				var workflowList = currentTarget.parentNode;
			}
		);

		return dynamicMenu;
	},



//// FETCH workflowTab JSON
//var splitLayout = '';
//dojo.xhrGet(
//	{
//		url: "../plugins/workflow/json/workflowTab.json",
//		handleAs: "json-comment-optional",
//		sync: true,
//		handle: function(response){
//			splitLayout = response.data;
//		}
//	}
//);


	splitLayout : {
		widgetType: "TabContainer",
		params: {id: "rootWidget", title: "Workflow 1", style:'height:450px; width:850px;'},
		children:
		[
			{
				widgetType: "SplitContainer",
				params: { title: 'Workflow1', id: "SplitContainer", orientation: "horizontal"},
				style: "width: 95%; height: 400px;",
				//style: "border: 3px solid grey; width: 95%; height: 400px;",
				children:
				[
					{
						widgetType: "ContentPane",
						params: {id: "workflow1leftPane", sizeShare:  15, style: "background: #EEEEFF;" },
						style: "background: #EEEEFF;",
						innerHTML: ""
					},
					{
						widgetType: "ContentPane",
						params: {id: "workflow1middlePane", sizeShare:  20, style: "background: #EEEEFF;" },
						style: "background: #EEEEFF;",
						innerHTML: ""
					},
					{
						widgetType: "ContentPane",
						params: {id: "workflow1rightPane", sizeShare:  75 },
						style: "background: #EEEEFF;",
						innerHTML: ""		
					}
				]
			}
		]
	},



	// GET APPLICATION INFO AND DISPLAY IN GRID IN RIGHT PANE
	loadApplicationInfo : function (event)
	{		

		// GET applicationName FROM NODE innerHTML
		var applicationName = event.target.innerHTML;

		// GET infoId FROM NODE 'infoId' ATTRIBUTE
		var infoId = event.target.getAttribute('infoId');

		// CONVERT APPLICATION NAME TO LOWER CASE
		var lowercaseApplicationName = applicationName.toLowerCase();

		// GET APPLICATION JSON
		var applicationData = '';	
		dojo.xhrGet(
			{
				url: "../plugins/workflow/json/application-" + lowercaseApplicationName + ".json",
				handleAs: "json",
					//handleAs: "json-comment-optional",
				sync: true,
				handle: function(response){
					applicationData = response.data;
				}
			}
		);

		// STORE IN WORKFLOW APPLICATIONS
		if ( ! this.workflowApplications[lowercaseApplicationName] )
		{

		}


		// SET APPLICATION MODEL USING APPLICATION DATA
		var applicationModel = new dojox.grid.data.Table(null, applicationData);

		// SET APPLICATION LAYOUT
		var applicationLayout = [
			{
				type: 'dojox.GridRowView', width: '0px'
			},
			{
				defaultCell: { width: 8, styles: 'text-align: center;'  },
				rows: [
					[
						{ name: 'Item', width: 15 },
						{ name: 'Info', editor: dojox.grid.editors.Input, styles: '', width: 20 }
					]
				]
			}
		];

		// GET GRID CONTAINER
		var gridContainer = dojo.byId(infoId);

		// REMOVE EXISTING GRID NODE FROM GRID CONTAINER
		while(gridContainer.firstChild)
		{
		   gridContainer.removeChild(gridContainer.firstChild);
		}

		// GET NEW UNIQUE GRID ID	
		var gridId = dojo.dnd.getUniqueId();

		// CREATE GRID
		var grid = new dojox.Grid({
			"id": gridId,
			"model": applicationModel,
			"structure": applicationLayout
		});

		// APPEND GRID TO GRID CONTAINER
		gridContainer.appendChild(grid.domNode);

		// NB: MUST ADD CLASS 'work-pane-split' TO SPLIT CONTAINER
		// IN ORDER FOR IT TO BE VISIBLE
		dojo.addClass(gridContainer, 'work-pane-split');

		// RENDER GRID
		grid.render();
	},


	// LOAD THE LEFT, MIDDLE AND RIGHT PANES OF THE WORKFLOW TAB
	loadWorkflowTab : function ()
	{

		var sourceId = this.paneId + "leftPane";
		var targetId = this.paneId + "middlePane";
		var infoId = this.paneId + "rightPane";


		var dragSource = new dojo.dnd.Source(
			sourceId,
			{
				accept: ["none"],
				copyOnly: true,
			}
		);

		var applications = [
			{ data: "Base callers",     type: ["notDraggable"]    },
			{ data: "Bustard",          type: ["isDraggable"]    },
			{ data: "Alta-Cyclic",      type: ["isDraggable"]    },
			{ data: "Assemblers",       type: ["notDraggable"]    },
			{ data: "AMOS",             type: ["isDraggable"]    },
			{ data: "Eland",            type: ["isDraggable"]    },
			{ data: "MAQ",              type: ["isDraggable"]    },
			{ data: "SOAP",             type: ["isDraggable"]    },
			{ data: "Velvet",           type: ["isDraggable"]    }
		];
		dragSource.insertNodes(false, applications);

		// GET this FOR USE IN LOOP BELOW
		var workflow = this;

		var allNodes = dragSource.getAllNodes();
		for ( var i = 0; i < allNodes.length; i++ )
		{
			var node = allNodes[i];
			var nodeClass = applications[i].type;
			var applicationName = node.innerHTML;

			// ADD CLASS FROM type TO NODE
			dojo.addClass(node, nodeClass);

			// ADD infoId TO NODE
			node.setAttribute('infoId', infoId);


			// IF ITEM CLASS IS isDraggable, SHOW APPLICATION INFO WHEN CLICKED
			if ( nodeClass == 'isDraggable' )
			{

				dojo.connect(node, "onclick",  null, function(e)
					{
						workflow.loadApplicationInfo(e, infoId);
					}
				);
			}

			//// IF ITEM isDraggable, SHOW APPLICATION INFO WHEN CLICKED
			//if ( nodeClass == 'isDraggable' )
			//{
			//
			//	allNodes[i].onclick = function()
			//	{
			//		workflow.loadApplicationInfo(applicationName, infoId);
			//	}
			//}

		}


		// SET NAMES FOR dropTargetNode AND ITS CONTAINER
		var dropTargetId = targetId + "-listNode";
		var dropTargetContainerId = targetId + "-container"; 	

		// MAKE CONTAINER FOR dropTarget ORDERED LIST
		var dropTargetContainer = document.createElement('div');
		dropTargetContainer.id = dropTargetContainerId;

		// CREATE targetNode ORDERED LIST NODE
		var dropTargetTitle = document.createElement('div');

		var dropTargetNode = document.createElement('ol');
		dojo.addClass(dropTargetNode, "dropTarget");
		dropTargetNode.innerHTML = "<span class='orderedList'> Workflow </span>";
		dropTargetNode.id = dropTargetId;

		// APPEND ORDERED LIST NODE TO CONTAINER
		dropTargetContainer.appendChild(dropTargetNode);

		// CREATE DND dropTarget OBJECT		
		var dropTarget = new dojo.dnd.Target( dropTargetNode, { accept: ["isDraggable"] } );

		// SAVE this AS WORKFLOW FOR USE INSIDE OVERRIDE OF onDndDrop BELOW
		var workflow = this;

		// OVERRIDE onDndDrop TO USE dojo.connect TO ADD EVENT TO NEW ITEM
		dropTarget.onDndDrop = function(source, nodes, copy)
		{
			// summary: topic event processor for /dnd/drop, called to finish the DnD operation
			// source: Object: the source which provides items
			// nodes: Array: the list of transferred items
			// copy: Boolean: copy items, if true, move items otherwise

			//break box
			do
			{ 

				if ( this.containerState != "Over" )
				{
					//this.onDndCancel();
					break;
				}
				var oldCreator = this._normalizedCreator;
				if(this != source){

					// transferring nodes from the source to the target

					if (this.creator)
					{
						// use defined creator
						this._normalizedCreator = function(node, hint){
							return oldCreator.call(this, source.getItem(node.id).data, hint);
						};
					}
					else
					{
						// we have no creator defined => move/clone nodes
						if(copy)
						{
							// clone nodes
							this._normalizedCreator = function(node, hint)
							{
								var t = source.getItem(node.id);
								var n = node.cloneNode(true);
								n.id = dojo.dnd.getUniqueId();
								var nodeId = n.id;

								// ADD ONCLICK TO LOAD APPLICATION INFO
								n.onclick = function(e)
								{
									if ( n.className.match(/isDraggable/) )
									{
										//Workflow.loadApplicationInfo(n.innerHTML, infoId);
										workflow.loadApplicationInfo(e, infoId);
									}

								}

								return {node: n, data: t.data, type: t.type};
							};

						}
						else
						{


							// move nodes
							this._normalizedCreator = function(node, hint){

								var t = source.getItem(node.id);
								source.delItem(node.id);
								return {node: node, data: t.data, type: t.type};
							};
						}
					}
				}
				else
				{

					// transferring nodes within the single source
					if(this.current && this.current.id in this.selection){ break; }

					if(this.creator)
					{
						// use defined creator
						if (copy)
						{
							// create new copies of data items
							this._normalizedCreator = function(node, hint){
								return oldCreator.call(this, source.getItem(node.id).data, hint);
							};
						}else
						{
							// move nodes
							if(!this.current){ break; }
							this._normalizedCreator = function(node, hint){
								var t = source.getItem(node.id);
								return {node: node, data: t.data, type: t.type};
							};
						}
					}
					else
					{
						// we have no creator defined => move/clone nodes


						if(copy)
						{
							// clone nodes
							this._normalizedCreator = function(node, hint){
								var t = source.getItem(node.id);
								var n = node.cloneNode(true);
								n.id = dojo.dnd.getUniqueId();
								return {node: n, data: t.data, type: t.type};
							};
						}
						else
						{
							// move nodes
							if(!this.current){ break; }
							this._normalizedCreator = function(node, hint){
								var t = source.getItem(node.id);
								return {node: node, data: t.data, type: t.type};
							};
						}
					}
				}
				this._removeSelection();
				if(this != source){
					this._removeAnchor();
				}
				if(this != source && !copy && !this.creator){
					source.selectNone();
				}
				this.insertNodes(true, nodes, this.before, this.current);
				if(this != source && !copy && this.creator){
					source.deleteSelectedNodes();
				}
				this._normalizedCreator = oldCreator;
			}
			while(false);
			// end of 'do'

			// GET JSON FOR APPLICATION INFO FOR NEWLY ADDED NODE
			// BY RUNNING THROUGH ALL THE EXISTING NODES AND COMPARING
			// WITH PREVIOUSLY EXISTING NODE LIST


			var elementArray = dojo.byId(workflow.paneId + 'middlePane').getElementsByTagName('div') ;
			for ( var i = 0; i < elementArray.length ; i++ )
			{
			}



			this.onDndCancel();
		}		


		// APPEND TARGET CONTAINER TO TARGET 
		dojo.byId(targetId).appendChild(dropTargetContainer);


		// GENERATE DYNAMIC DHTML MENU
		var dynamicMenu = this.createMenu();

		// BIND THE MENU TO THE DND NODES SO WE OPEN THE
		// MENU WHEN WE CLICK OVER THE NODES
		dynamicMenu.bindDomNode( dojo.byId(sourceId) );
		dynamicMenu.bindDomNode( dojo.byId(targetId) );


	} // end of Workflow.loadWorkflowTab


}); // end of Workflow	

