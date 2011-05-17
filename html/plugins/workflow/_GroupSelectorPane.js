
// HANDLE GETTING SELECTED OBJECT IN DND LIST
// OVERRIDE onDropInternal(), onDropExternal()
// OR onDrop TO COVER BOTH ON THE source/target OBJECT.

// TEST URL:
//http://localhost:8080/Bioptic0.2.3/html/file-drag2.html

dojo.provide("plugins.workflow._GroupSelectorPane");

// DEPENDENCIES
dojo.require("plugins.workflow._GroupDragPane");

dojo.declare("plugins.workflow._GroupSelectorPane",
	[plugins.workflow._GroupDragPane], {

	// callback FROM OBJECT THAT GENERATED THE FileManager WHICH CREATED
	// THE FileSelector THAT CREATED THIS _GroupSelectorPane
	callback : null,
	parentWidget : null,
	grandParentWidget : null,
	workflowObject : null,

	// SET callingObject IN CONSTRUCTOR
	constructor : function(args) {


		this.callback = args.callback;
		this.parentWidget = args.parentWidget;
		this.grandParentWidget = args.grandParentWidget;
		this.workflowObject = args.workflowObject;



	//	this.inherited(arguments);
	},



	// ADD PROGRAMMATIC CONTEXT MENU
	createMenu : function (type)
	{

		var selectorPane = this;

		var node = document.createElement('div');
		var nodeId = dojo.dnd.getUniqueId();
		node.id = nodeId;
		var dynamicMenu = new dijit.Menu( { id: nodeId	} );

		// STUB FOR SETTING WORKFLOW-SPECIFIC MENUS
		if ( type == "workflow" )
		{

		}

		// ADD MENU TITLE
		dynamicMenu.addChild(new dijit.MenuItem( { label:"Application Menu", disabled: false } ));
		dynamicMenu.addChild(new dijit.MenuSeparator());

		//// ONE OF FOUR WAYS TO DO MENU CALLBACK WITH ACCESS TO THE MENU ITEM AND THE CURRENT TARGET 	
		// 4. dojo.connect CALL
		//	REQUIRES:
		//		ADDED menu.currentTarget SLOT TO dijit.menu
		var mItem1 = new dijit.MenuItem(
			{
				//id: "workflow" + this.paneId + "remove",
				label: "Select",
				disabled: false
			}
		);
		dynamicMenu.addChild(mItem1);


		dojo.connect(mItem1, "onClick", function()
			{

				// GET FULL PATH OF FILE
				var item = dynamicMenu.currentTarget.item;
				var file = item.parentPath + "/" + item.path;


				// USE LOCATION IN THE CASE OF SHARED FILES
				var location = selectorPane.parentWidget.location;

				// SET TYPE (FILE OR DIRECTORY)
				var type = "file";
				if ( item.directory == true )	type = "directory";

				// DO CALLBACK TO CHANGE OPTION VALUE
				selectorPane.callback(file, location, type);

				// SET 'file' IN CALL ARGUMENTS
				//selectorPane.callback.callArguments.file = file;
				//selectorPane.callback.callArguments.location = location;

				// SET grandParentWidget IN CALL ARGUMENTS
				//selectorPane.callback.callArguments.grandParentWidget = selectorPane.grandParentWidget;
				//selectorPane.callback.callArguments.workflowObject = selectorPane.workflowObject;

				// DESTROY GRAND PARENT WIDGET (FILE MANAGER)
				selectorPane.grandParentWidget.fileManagerDialog.hide();

				// DO CALLBACK	
			//	selectorPane.callback.callFunction(selectorPane.callback.callArguments);

			}
		);

		return dynamicMenu;
	}



});

