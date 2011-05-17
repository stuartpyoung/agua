
// HANDLE GETTING SELECTED OBJECT IN DND LIST
// OVERRIDE onDropInternal(), onDropExternal()
// OR onDrop TO COVER BOTH ON THE source/target OBJECT.

// TEST URL:
//http://localhost:8080/Bioptic0.2.3/html/file-drag2.html

dojo.provide("plugins.files._GroupDragPane");

// REQUIRE EXTERNAL MODULES
// Drag n Drop
dojo.require("dojo.dnd.Source");
// Menu
dojo.require("dijit.Menu");

// FILE UPLOAD
dojo.require("plugins.files.FileUpload");

// POPUP DIALOGUE
dojo.require("plugins.files.Dialog");

// RENAME AND NEW FOLDER
dojo.require("dijit.InlineEditBox");

dojo.declare("plugins.files._GroupDragPane",
	[dojox.widget._RollingListPane], {

	// summary: a pane that will handle groups (treats them as menu items)

	// templateString: string
	//	our template
	templateString: '<div><div dojoAttachPoint="containerNode"></div>' +
					'<div dojoAttachPoint="menuContainer">' +
						'<div dojoAttachPoint="menuNode"></div>' +
					'</div></div>',

	// _dragSource: dijit.Menu
	//  The menu that we will call addChild() on for adding items
	_dragSource: null,

	// dialog box to show copying
	_copyBox : null,

	// OBJECT-WIDE DEBUG STATUS
	debug : false,


	_doQuery: function(){


		// summary: either runs the query or loads potentially not-yet-loaded items.
		this.isLoaded = false;

			this._setContentAndScroll(this.onFetchStart());


			this.store.fetch(
				{
					query: this.store.path, 
					onComplete: function(items){
						this.items = items;
						this.onItems();
					}, 
					onError: function(e){
						this._onError("Fetch", e);
					},
					scope: this
				}
			);
		//}
	},



_onLoadHandler: function(data){
// OVERRIDE TO AVOID THIS ERROR:
// Error undefined running custom onLoad code: This deferred has already been resolved


if ( data == null )
{
}


	// summary:
	//		This is called whenever new content is being loaded
	this.isLoaded = true;
	try{
		//this.onLoadDeferred.callback(data);
		this.onLoad(data);
	}catch(e){
		console.error('Error '+this.widgetId+' running custom onLoad code: ' + e.message);
	}
},




_loadCheck: function(/* Boolean? */ forceLoad){


	// summary: checks that the store is loaded
	var displayState = this._isShown();
	if((this.store || this.items) && (forceLoad || (this.refreshOnShow && displayState) || (!this.isLoaded && displayState))){

		this.query = this.path;
		//var request={};
		//request.query = this.path;
		//this._doQuery();

		this._doQuery();
			//{
			//	query: this.path
			//}
		//);

	}
},

_setContent: function(/*String|DomNode|Nodelist*/cont){
	if(!this._dragSource){
		// Only set the content if we don't already have a menu
		this.inherited(arguments);
	}
},

postCreate: function(){
//	Called after a widget's dom has been setup
},

createMenu : function (type) {
// ADD PROGRAMMATIC CONTEXT MENU


	if ( type == "workflow" )
	{
		if ( this.workflowMenu != null ) return this.workflowMenu.menu;
	}
	else if ( type == "folder" )
	{
		if ( this.folderMenu != null ) return this.folderMenu.menu;
	}
	else
	{
		if ( this.fileMenu != null ) return this.fileMenu.menu;
	}
},


onItems : function() {
//	called after a fetch or load - at this point, this.items should be
//  set and loaded.

	var groupDragPane = this;

	for ( var i = 0; i < this.items.length; i++)
	{			
		this.items[i].path = this.items[i].name;
	}	

	var selectItem, hadChildren = false;

	this._dragSource = this._getDragSource(
		{
			path: this.path,
			parentPath: this.parentPath
		}
	);

	// ADD THE STORE'S parentPath TO THE MENU
	this._dragSource.store = this.store;

	// IF THERE ARE NO ITEMS FOR THIS DIRECTORY,
	// INSERT A FAKE ITEM WITH NAME <EMPTY> 
	if ( ! this.items.length )
	{

		// THIS CHILD SHOULD READ "<EMPTY>"
		var item = new Object;
		item.name = "<EMPTY>";
		item.parentPath = this.parentPath;
		this.items.push(item);
	}

	var child, selectMenuItem;
	if ( this.items.length )
	{
		dojo.forEach(
			this.items,
			function(item)
		{
			child = this.parentWidget._getMenuItemForItem(item, this);

			if ( child.data == 'EMPTY' )
			{
			}

			if ( child )
			{
				this._dragSource.insertNodes(false, [ child ]);

				var insertedNodes = this._dragSource.getAllNodes();
				var lastNode = insertedNodes[insertedNodes.length - 1];

				// ADD DATA STORE TO ITEM'S CHILDREN
				if ( item.children )
				{
					for ( var i = 0; i < item.children.length; i++ )
					{
						item.children[i]._S = item._S;
						var childParentPath;

						var fullPath = '';
						if ( item.parentPath )
						{
							fullPath = item.parentPath;
						}

						if ( item.path )
						{
							fullPath += "/" + item.path;

						}

						if ( fullPath )
						{
							item.children[i].parentPath = fullPath;
						}
					}
					item.directory = true;
				}

				// ADD ITEM TO THIS NODE
				lastNode.item = item;

				// SET CLASS
				dojo.hasClass(lastNode, "dojoxRollingListItemSelected");

				var applicationName = lastNode.innerHTML;

				// GET indexInParent - THE LEVEL OF THIS MENU IN THE PARENT
				var indexInParent = this.getIndexInParent();

				// SET nodeType BASED ON THE indexInParent TO COINCIDE WITH accept PARAMETER
				// OF DND SOURCE OF SAME LEVEL (E.G., Workflows CAN BE DRAGGED NEXT TO OTHER
				// WORKFLOWS BUT NOT INTO THE LOWER FILE DIRECTORIES)
				var nodeType;
				if ( indexInParent == 0 )
				{
					nodeType = 'workflow';
				}
				else
				{
					if ( item.directory == false )
						nodeType = "file";
					else
						nodeType = "folder";
				}

				// GENERATE DYNAMIC DHTML MENU
				var dynamicMenu = groupDragPane.createMenu(nodeType);

				// BIND THE MENU TO THE DND NODE
				if ( dynamicMenu != null )
				{
					dynamicMenu.bindDomNode( lastNode );
				}
				lastNode.setAttribute("dndType", nodeType);					

				// CONNECT ONCLICK
				dojo.connect(lastNode, "onclick", this, "onclickHandler");

				// OVERRIDE dojo.dnd.Source.onDropExternal TO NOTIFY SERVER OF CHANGES

				//this._dragSource.onDropExternal = function(e) {
				this._dragSource.onDropExternal =  function(source, nodes, copy){						


					var dragSource = this;

/*
					// COMPLETE THE COPY ON THE REMOTE FILESYSTEM
					//
					// 1. CARRY OUT DND COPY
					//
					// 2. CHECK IF FILENAME ALREADY EXISTS, IF SO
					// 		DO POPUP TO CONFIRM OVERWRITE	
					//
					// 3. MESSAGE SERVER TO COPY FILES
					//
					// 4. SHOW ANIMATED 'COPYING' DIALOGUE
					//
					// 5. POLL SERVER FOR STATUS AND WAIT UNTIL COMPLETE
					//
					// 6. RELOAD THE PANE TO SHOW THE NEW FILE SYSTEM
*/


/*					// 2. CHECK IF FILENAME ALREADY EXISTS
					//
					//		IF SO, DO POPUP TO CONFIRM OVERWRITE
					//
					//		OTHERWISE, MESSAGE SERVER WITH COPY
					//
*/
					var copyFile = function (url, query)
					{
						dojo.xhrGet(
							{
								url: url + query,
								handleAs: "json",
									//handleAs: "json-comment-optional",
								sync: false,
								handle: function(response){
									if ( response.error )
									{
										fileCopyError = response.error;

										// POP UP OVERWRITE CONFIRMATION DIALOGUE WINDOW
										if ( response.error.match(/^File exists$/) )
										{
											checkOverwrite(url, query);
										}
										else
										{
											return;
										}
									}
									if ( ! response.error )
									{
										doCopy(url, query);
									}

								}
							}
						);
					}


					//// EXTRACT owner FROM STORE url
					var url = nodes[0].item._S.url;
					var owner;
					if ( url.match(/owner=([^&]+)/) )
					{
						owner = url.match(/owner=([^&]+)/)[1];
					}
					var groupname;
					if ( url.match(/groupname=([^&]+)/) )
					{
						groupname = url.match(/groupname=([^&]+)/)[1];
					}



					// 3. MESSAGE SERVER TO COPY FILES
					var url = Agua.cgiUrl + "project.cgi?";
					var file = nodes[0].item.parentPath + "/" + nodes[0].item.path;
					var destination = this.path;


					//var query = "mode=copyFile&sessionId=1228791394.7868.158&username=admin&file=" + file + "&destination=" + destination;
					var query = "mode=copyFile";
					query += "&sessionId=" + Agua.cookie('sessionId');
					query += "&username=" + Agua.cookie('username');
					query += "&file=" + file;
					query += "&destination=" + destination;
					if ( owner != null && owner != '' )
					{
						query += "&owner=" + owner;

					}
					if ( groupname != null && groupname != '' )
					{
						query += "&groupname=" + groupname;

					}
					copyFile(url, query);

					var fileCopyError;

					var itemParent = function(item)
					{
						// SET DIRECTORY = TRUE
						item.directory = true;

						// CHANGE NAME, PATH AND PARENTPATH
						//
						// 1. IF PARENTPATH CONTAINS MULTIPLE LEVELS, E.G., 'Project1/Workflow1-assembly',
						if ( item.parentPath.match(/^.+\/([^\/]+)$/) )
						{
							item.path = item.parentPath.match(/^(.+?)\/([^\/]+)$/)[2];


							item.parentPath = item.parentPath.match(/^(.+?)\/([^\/]+)$/)[1];
						}

						// 2. IF PARENTPATH IS AT THE TOP LEVEL, E.G., 'Project1',
						// 		SET PATH = PARENTPATH AND PARENTPATH = '' 
						else if ( item.parentPath.match(/\/*[^\/]+$/) )
						{
							item.path = item.parentPath;
							item.parentPath = '';
						}

						item.name = item.path;

						//item.parentPath = item.parentPath.match(/^(.+?)\/[^\/]+$/)[1];

						return item;
					}

					var reloadPane = function()
					{

						var item = lastNode.item;

						var children = item.children;
						if ( ! children )
						{
							children = item.items;
						}

						// CHANGE item PATH, NAME AND PARENTPATH TO ONE FOLDER UP
						item = itemParent(item);

						var itemPane = groupDragPane.parentWidget._getPaneForItem(item, this, children);

						groupDragPane.query = itemPane.store.query;

						if(itemPane)
						{
							// CALLS addChild IN plugins.files._GroupDragPane
							// summary: adds a child to this rolling list - if passed an insertIndex,
							//  then all children from that index on will be removed and destroyed
							//  before adding the child.

							if ( itemPane.store && itemPane.store.path )
							{
							}

							var paneIndex = groupDragPane.getIndexInParent();

							//groupDragPane.parentWidget.addChild( itemPane, groupDragPane.getIndexInParent() + 1);
							groupDragPane.parentWidget.addChild( itemPane, groupDragPane.getIndexInParent() );
						}
						//else
						//{
						//	groupDragPane.parentWidget(groupDragPane);
						//
						//	groupDragPane.parentWidget._onItemClick(null, groupDragPane, selectMenuItem.item, selectMenuItem.children);
						//}

					};


					var checkOverwrite = function(url, query)
					{

						var fileDragTop = groupDragPane.parentWidget.domNode.offsetTop;
						var fileDragLeft = groupDragPane.parentWidget.domNode.offsetLeft;
						var fileDragHeight = groupDragPane.parentWidget.domNode.offsetHeight;
						var fileDragWidth = groupDragPane.parentWidget.domNode.offsetWidth;

						var dialogWindow = new plugins.files.Dialog( { url: url, query: query, parentFunction: copyFile, message: "file " + file + " already exists in folder: " + destination + "\nDo you want to overwrite it?"} );

					}; // var checkOverwrite = function(url, query)


					// 4. SHOW ANIMATED 'COPYING' DIALOGUE
					var doCopy = function(url, query)
					{
						dojo.require("dojo._base.fx");
						dojo.require("dojo._base.html");


						var fileDragTop = groupDragPane.parentWidget.domNode.offsetTop;
						var fileDragLeft = groupDragPane.parentWidget.domNode.offsetLeft;
						var fileDragHeight = groupDragPane.parentWidget.domNode.offsetHeight;
						var fileDragWidth = groupDragPane.parentWidget.domNode.offsetWidth;

						// CONSTRUCT A DIALOGUE BOX THE SAME SIZE AS THE PANE						
						var box = document.createElement('div');
						groupDragPane.parentWidget.domNode.appendChild(box);

						// SAVE BOX IN GROUPDRAGPANE
						groupDragPane._copyBox = box;

						//document.body.appendChild(box);
						box.setAttribute('class', 'copyDialog');
						//box.style.top =  fileDragTop + "px";
						//box.style.left = fileDragLeft + "px";
						dojo.style(
							box,
							{
								"opacity": "0.3",
								"position": "absolute",
								"top": fileDragTop + "px",
								"left": fileDragLeft + "px"
							}
						);

						// CREATE DIV FOR TITLE
						var titleDiv = document.createElement('div');
						box.appendChild(titleDiv);
						titleDiv.setAttribute('class', 'titleDiv');

						// CREATE TITLE TEXT
						var title = document.createTextNode("Copying...");
						titleDiv.appendChild(title);

						var counter = 0;
						var copyCompleted = false;

						var handleCopy = function(response)
						{
							if ( response.status == 'completed' )
							{
								copyCompleted = true;
							}
						};

						var copyDialog = function (url, query, callbackFunction)
						{

							dojo.style(
								box,
								{
									"opacity": "0.3"
								}
							);

							dojo.fadeIn(
								{
									node: box,
									delay: 1,
									duration: 1000,
									properties: {
										opacity: {
											start: 1.0,
											end: 0.0
										},
										width:
										{
											start: fileDragWidth,
											end: fileDragWidth
										},
										height:
										{
											start: fileDragHeight,
											end: fileDragHeight
										}
									},

									//easing: function(n){
									//	return (n==0) ? 0 : Math.pow(2, 10 * (n - 1));
									//},

									onEnd: function (){

										// 6. IF COPY IS COMPLETED, RELOAD THE PANE TO SHOW THE NEW FILE SYSTEM
										if ( copyCompleted )
										{
											var box = groupDragPane._copyBox;
											if ( box != null )
												groupDragPane.parentWidget.domNode.removeChild(box);

											//document.body.removeChild(box);
											reloadPane();

										}

										counter++;

										// 5. POLL SERVER FOR STATUS AND WAIT UNTIL COMPLETE
										if ( ! query.match(/modifier/) )
										{
											query += "&modifier=status";
										}
										else
										{
											query = query.replace(/modifier=overwrite/, "modifier=status");
										}


										dojo.xhrGet(
											{
												url: url + query,
												handleAs: "json",
													//handleAs: "json-comment-optional",
												sync: false,
												handle: callbackFunction
											}
										);

										if ( ! copyCompleted )
										{
											//// GIVE UP AFTER 8 GOES
											//if ( counter < 8 )
											//{
											//	copyDialog(url, query, callbackFunction);
											//}
										}
										else
										{
										}
									}
								}

							).play();

						}; // var copyDialog = function()

						copyDialog(url, query, handleCopy);

					} // var doCopy = function(url, query)



					// SAVE THIS - DECIDE LATER WHETHER TO KEEP 
					var fileDrop = function()
					{
						// 1. CARRY OUT DND COPY (FROM OVERRIDDEN FUNCTION)
						var oldCreator = dragSource._normalizedCreator;
						// transferring nodes from the source to the target
						if(dragSource.creator)
						{
							// use defined creator
							dragSource._normalizedCreator = function(node, hint)
							{
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
									return {node: n, data: t.data, type: t.type};
								};
							}
							else
							{
								// move nodes
								this._normalizedCreator = function(node, hint)
								{
									var t = source.getItem(node.id);
									source.delItem(node.id);
									return {node: node, data: t.data, type: t.type};
								};
							}
						}

						this.selectNone();
						if(!copy && !this.creator)
						{
							source.selectNone();
						}
						this.insertNodes(true, nodes, this.before, this.current);
						if(!copy && this.creator)
						{
							source.deleteSelectedNodes();
						}
						this._normalizedCreator = oldCreator;

					}

				} // this._dragSource.onDropExternal

			} // if(child)

		}, this); // dojo.forEach(this.items, function(item)

	} // if ( this.items.length )



	// ADD dojo.connect TO DND SOURCE NODES
	var allNodes = this._dragSource.getAllNodes();
	for ( var i = 0; i < allNodes.length; i++ )
	{
		var node = allNodes[i];
		dojo.addClass(node, "fileDrag");

		if ( node.item.directory )
		{
			dojo.addClass(node, "directory");
			//dojo.addClass(node, "directoryClosed");
		}
		else
		{
			dojo.addClass(node, "file");
			//dojo.addClass(node, "fileClosed");
		}

/*
		// GENERATE ICON NODE TO BE PLACED IN FRONT OF NODE
		var iconNode = document.createElement('div');
		node.appendChild(iconNode);

		//// ADD CLASS
		//dojo.removeClass(node, "dojoDndItem");
		//dojo.addClass(node, "soria dojoxFileDragItem dojoxDirectoryItemIcon");
		dojo.addClass(iconNode, "directoryIcon");
*/

		var applicationName = node.innerHTML;

		// GET indexInParent - THE LEVEL OF THIS MENU IN THE PARENT
		var indexInParent = this.getIndexInParent();

		// SET nodeType BASED ON THE indexInParent TO COINCIDE WITH accept PARAMETER
		// OF DND SOURCE OF SAME LEVEL (E.G., Workflows CAN BE DRAGGED NEXT TO OTHER
		// WORKFLOWS BUT NOT INTO THE LOWER FILE DIRECTORIES)
		var nodeType;
		if ( indexInParent == 0 )
		{
			nodeType = 'workflow';
		}
		else
		{
			nodeType = "file";
		}
		node.setAttribute("dndType", nodeType);
	}	

	// ADD CLASS - THIS ADDS IT TO THE CONTAINER OF THE DND SOURCE - NOT WHAT WE WANT
	//dojo.addClass(this.containerNode, "dojoxDirectoryItemIcon");
	//.dojoxDirectoryItemIcon


	this.containerNode.innerHTML = "";		
	this.containerNode.appendChild(this.menuNode);
	this.parentWidget.scrollIntoView(this);


return;

	this.inherited(arguments);

},

// HANDLE CLICK ON FILE OR FOLDER
onclickHandler : function(e) {

	// GET THE CLICKED DND SOURCE ITEM NODE
	var item = e.target.item;
	var children = item.children;
	if ( ! children )
	{
		children = item.items;
	}
	var itemPane = this.parentWidget._getPaneForItem(item, this, children);

	//
	//return;


	// SET this.query TO itemPane.store.query
	this.query = itemPane.store.query;

	if(itemPane)
	{
		// CALLS addChild IN plugins.files._GroupDragPane
		// summary: adds a child to this rolling list - if passed an insertIndex,
		//  then all children from that index on will be removed and destroyed
		//  before adding the child.

		if ( itemPane.store && itemPane.store.path )
		{
		}

		this.parentWidget.addChild(itemPane, this.getIndexInParent() + 1);
	}
	else
	{
		this.parentWidget(this);
		this.parentWidget._onItemClick(null, this, selectMenuItem.item, selectMenuItem.children);
	}
},



startup: function(){

	this.inherited(arguments);
	//this.parentWidget._updateClass(this.domNode, "GroupPane");
},



focus: function(/*boolean*/force){
	// summary: sets the focus to this current widget


	if(this._dragSource){
		if(this._pendingFocus){
			this.disconnect(this._pendingFocus);
		}
		delete this._pendingFocus;

		// We focus the right widget - either the focusedChild, the
		//   selected node, the first menu item, or the menu itself
		var focusWidget = this._dragSource.focusedChild;
		if(!focusWidget){
			var focusNode = dojo.query(".dojoxRollingListItemSelected", this.domNode)[0];
			if(focusNode){
				focusWidget = dijit.byNode(focusNode);
			}
		}

		if(!focusWidget){
			focusWidget = this._dragSource.getAllNodes()[0] || this._dragSource;
		}

		this._focusByNode = false;

		if(focusWidget.focusNode){
			if(!this.parentWidget._savedFocus || force){
				try{focusWidget.focusNode.focus();}catch(e){}
			}
			window.setTimeout(function(){
				try{
					dijit.scrollIntoView(focusWidget.focusNode);
				}catch(e){}
			}, 1);
		}else if(focusWidget.focus){
			if(!this.parentWidget._savedFocus || force){
				focusWidget.focus();
			}
		}else{
			this._focusByNode = true;
		}
		this.inherited(arguments);
	}else if(!this._pendingFocus){
		this._pendingFocus = this.connect(this, "onItems", "focus");
	}
	else
	{
	}

},





_getDragSource: function(){
	// summary: returns a widget to be used for the container widget.

	// GET UNIQUE ID FOR THIS MENU TO BE USED IN DND SOURCE LATER
	var objectName = "dojo.dnd.Source";
	var id = dijit.getUniqueId(objectName.replace(/\./g,"_"));
	//var id = dijit.getUniqueId(this.declaredClass.replace(/\./g,"_"));		

	// SET THE MENU NODE'S ID TO THIS NEW ID
	this.menuNode.id = id;

	// GET indexInParent - THE LEVEL OF THIS DRAG SOURCE IN THE PARENT
	var indexInParent = this.getIndexInParent();

	// SET accept BASED ON THE indexInParent

	var acceptType;
	if ( indexInParent == 0 )
	{
		acceptType = 'workflow';
	}
	else
	{
		acceptType = "file";
	}

	// GENERATE DND SOURCE WITH UNIQUE ID
	var dragSource = new dojo.dnd.Source(
		id,
		{
			accept: [ acceptType ],
			copyOnly: true
		}
	);


	// SET baseClass
	this.menuNode.setAttribute('class', 'fileDrag');

	// SET PARENTPATH AND PATH
	if ( this.path )
	{
		dragSource.path = this.path;
	}
	if ( this.parentPath )
	{
		dragSource.parentPath = this.parentPath;
	}


	if(!dragSource._started){
		dragSource.startup();
	}


	return dragSource;
},



// LEGACY FUNCTION FROM RollingListGroupPane. REVIEW WHETHER TO REMOVE LATER.
_getSelected: function(/*dijit.Menu?*/ menu){
	// summary:
	//	returns the selected menu item - or null if none are selected


	if( !menu )   { menu = this._dragSource; }

	if ( menu )
	{
		// SELECTED DND SOURCE HAS CLASS dojoDndItemAnchor
		var nodes = this._dragSource.getAllNodes();

//			var children = this._dragSource.getChildren();
		for ( var i = 0; i < nodes.length; i++ )
//			for(var i = 0, item; (item = children[i]); i++)
		{
			//if(dojo.hasClass(item.domNode, "dojoxRollingListItemSelected")){
			var node = nodes[i];
			if ( node.className.match(/dojoDndItemOver/) )
			{
				return node;
			}
		}
	}
	return null;
},



// LEGACY FUNCTION FROM RollingListGroupPane. REVIEW WHETHER TO REMOVE LATER.

_setSelected: function(/*dijit.MenuItem?*/ item, /*dijit.Menu?*/ menu){
	// summary:
	//	selects the given item in the given menu (defaults to pane's menu)
	if(!menu){ menu = this._dragSource;}
	if(menu){
		dojo.forEach(menu.getChildren(), function(i){
			this.parentWidget._updateClass(i.domNode, "Item", {"Selected": (item && (i == item && !i.disabled))});
		}, this);
	}
}
});
