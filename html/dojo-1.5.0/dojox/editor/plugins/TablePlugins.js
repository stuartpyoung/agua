dojo.provide("dojox.editor.plugins.TablePlugins");
dojo.require("dijit._editor._Plugin");
dojo.require("dijit._editor.selection");
dojo.require("dijit.Menu");
dojo.require("dojo.i18n");
dojo.requireLocalization("dojox.editor.plugins", "TableDialog");

dojo.experimental("dojox.editor.plugins.TablePlugins");

// summary:
//		A series of plugins that give the Editor the ability to create and edit 
//		HTML tables. See the end of this document for all avaiable plugins
//		and dojox/editorPlugins/tests/editorTablePlugs.html for an example
//
// example:
//		|	<div dojoType="dijit.Editor" plugins="[
//		|			'bold','italic','|',
//		|			{name: 'dojox.editor.plugins.TablePlugins', command: 'insertTable'},
//		|			{name: 'dojox.editor.plugins.TablePlugins', command: 'modifyTable'}
//		|		]">
//		|		Editor text is here
//		|	</div>
//
// TODO:	
//		Currently not supporting merging or splitting cells
//
// FIXME:	Undo is very buggy, and therefore unimeplented in all browsers 
//			except IE - which itself has only been lightly tested.
//
// FIXME:	Selecting multiple table cells in Firefox looks to be impossible.
//			This affect the 'colorTableCell' plugin. Cells can still be 
//			colored individually or in rows.

dojo.declare("dojox.editor.plugins._TableHandler", dijit._editor._Plugin,{
	// summary:
	//		A global object that handles common tasks for all the plugins. Since 
	//		there are several plugins that are all calling common methods, it's preferable
	//		that they call a centralized location that either has a set variable or a 
	//		timeout to only repeat code-heavy calls when necessary.
	//
	tablesConnected:false,
	currentlyAvailable: false,
	alwaysAvailable:false,
	availableCurrentlySet:false,
	initialized:false,
	tableData: null,
	shiftKeyDown:false,
	editorDomNode: null,
	undoEnabled: true, //Using custom undo for all browsers.
	refCount: 0, 

	doMixins: function(){

		dojo.mixin(this.editor,{
			getAncestorElement: function(tagName){
				return dojo.withGlobal(this.window, "getAncestorElement",dijit._editor.selection, [tagName]);
			},
			hasAncestorElement: function(tagName){
				return dojo.withGlobal(this.window, "hasAncestorElement",dijit._editor.selection, [tagName]);
			},
			selectElement: function(elem){
				dojo.withGlobal(this.window, "selectElement",dijit._editor.selection, [elem]);
			},
			byId: function(id){
				return dojo.withGlobal(this.window, "byId", dojo, [id]);
			},
			query: function(arg, scope, returnFirstOnly){
				// this shortcut is dubious - not sure scoping is necessary
				var ar = dojo.withGlobal(this.window, "query", dojo, [arg, scope]);
				return (returnFirstOnly) ? ar[0] : ar;
			}
		});

	},
	initialize: function(editor){
		// summary:
		//		Initialize the global handler upon a plugin's first instance of setEditor
		//

		// All plugins will attempt initialization. We only need to do so once.
		// But keep track so that it is cleaned up when all usage of it for an editor has
		// been removed.
		this.refCount++;

		// Turn on custom undo for all.
		editor.customUndo = true;

		if(this.initialized){ return; }

		this.initialized = true;
		this.editor = editor;

		this.editor._tablePluginHandler = this;

		//Editor loads async, can't assume doc is ready yet.  So, use the deferred of the
		//editor to init at the right time.
		editor.onLoadDeferred.addCallback(dojo.hitch(this, function(){
			this.editorDomNode = this.editor.editNode || this.editor.iframe.document.body.firstChild;

			// RichText should have a mouseup connection to recognize drag-selections
			// Example would be selecting multiple table cells
			this._myListeners = [];
			this._myListeners.push(dojo.connect(this.editorDomNode , "mouseup", this.editor, "onClick")); 
			this._myListeners.push(dojo.connect(this.editor, "onDisplayChanged", this, "checkAvailable"));
			this._myListeners.push(dojo.connect(this.editor, "onBlur", this, "checkAvailable"));
			this.doMixins();
			this.connectDraggable();
		}));
	},

	getTableInfo: function(forceNewData){
		// summary:
		//		Gets the table in focus
		//		Collects info on the table - see return params
		//
		if(forceNewData){ this._tempStoreTableData(false); }
		if(this.tableData){
			// tableData is set for a short amount of time, so that all 
			// plugins get the same return without doing the method over
			return this.tableData;	
		}
		var tr, trs, td, tds, tbl, cols, tdIndex, trIndex;

		td = this.editor.getAncestorElement("td");
		if(td){ tr = td.parentNode; }

		tbl = this.editor.getAncestorElement("table");

		tds = dojo.query("td", tbl);
		tds.forEach(function(d, i){
			if(td==d){tdIndex = i;}
		});
		trs = dojo.query("tr", tbl);
		trs.forEach(function(r, i){
			if(tr==r){trIndex = i;}
		});
		cols = tds.length/trs.length;
		var o = {
			tbl:tbl,		// focused table
			td:td,			// focused TD
			tr:tr,			// focused TR
			trs:trs,		// rows
			tds:tds,		// cells
			rows:trs.length,// row amount
			cols:cols,		// column amount
			tdIndex:tdIndex,// index of focused cell
			trIndex:trIndex,	// index of focused row
			colIndex:tdIndex%cols
		};
		this.tableData = o;
		this._tempStoreTableData(500);	
		return this.tableData;
	},

	connectDraggable: function(){
		// summary:
		//		Detects drag-n-drop in the editor (could probably be moved to there)
		//		Currently only checks if item dragged was a TABLE, and removes its align attr
		//		DOES NOT WORK IN FF - it could - but FF's drag detection is a monster
		//
		if(!dojo.isIE){
			//console.warn("Drag and Drop is currently only detectable in IE.");
			return;
		}

		// IE ONLY
		this.editorDomNode.ondragstart = dojo.hitch(this, "onDragStart");
		this.editorDomNode.ondragend = dojo.hitch(this, "onDragEnd");

		//NOTES:
		// FF _ Able to detect the drag-over object (the editor.domNode)
		//	Not able to detect an item's ondrag() event
		//	Don't know why - I actually got it working when there was an error
		//	Something to do with different documents or windows I'm sure
		//
		/*tbl.ondragstart=dojo.hitch(this, "onDragStart");

		tbl.addEventListener("dragstart", dojo.hitch(this, "onDragStart"), false);
		tbl.addEventListener("drag", dojo.hitch(this, "onDragStart2"), false);
		tbl.addEventListener("dragend", dojo.hitch(this, "onDragStart3"), false);

		dojo.withGlobal(this.editor.window, "selectElement",dijit._editor.selection, [tbl]);

		tbl.ondragstart = function(){
		};
		tbl.ondrag = function(){
			alert("drag")
		*/		
	},
	onDragStart: function(){
		var e = window.event;
		if(!e.srcElement.id){
			e.srcElement.id = "tbl_"+(new Date().getTime());	
		}
	},
	onDragEnd: function(){
		// summary:
		//		Detects that an object has been dragged into place
		//		Currently, this code is only used for when a table is dragged
		//		and clears the "align" attribute, so that the table will look
		//		to be more in the place that the user expected.
		//		TODO: This code can be used for other things, most 
		//		notably UNDO, which currently is not quite usable.
		//		This code could also find itself in the Editor code when it is 
		//		complete.

		var e = window.event;
		var node = e.srcElement;
		var id = node.id;
		var win = this.editor.window;

		// clearing a table's align attr
		// TODO: when ondrag becomes more robust, this code block
		//	should move to its own method
		if(node.tagName.toLowerCase()=="table"){
			setTimeout(function(){
				var node =  dojo.withGlobal(win, "byId", dojo, [id]);
				dojo.removeAttr(node, "align");
			}, 100);
		}
	},
	checkAvailable: function(){
		// summary:
		//		For table plugs
		//		Checking if a table or part of a table has focus so that 
		//		Plugs can change their status
		//
		if(this.availableCurrentlySet){
			// availableCurrentlySet is set for a short amount of time, so that all 
			// plugins get the same return without doing the method over
			return this.currentlyAvailable;
		}

		if(!this.editor) {
			return false;
		}
		if(this.alwaysAvailable) {
			return true;
		}

		// Only return avalable if the editor is focused.
		this.currentlyAvailable = this.editor._focused ? this.editor.hasAncestorElement("table") : false;

		if(this.currentlyAvailable){
			this.connectTableKeys();
		}else{
			this.disconnectTableKeys();
		}

		this._tempAvailability(500);
		dojo.publish(this.editor.id + "_tablePlugins", [ this.currentlyAvailable ]);
		return this.currentlyAvailable;
	},

	_prepareTable: function(tbl){
		//	For IE's sake, we are adding IDs to the TDs if none is there
		//	We go ahead and use it for other code for convenience
		//	
		var tds = this.editor.query("td", tbl);
		if(!tds[0].id){
			tds.forEach(function(td, i){
				if(!td.id){
					td.id = "tdid"+i+this.getTimeStamp();
				}
			}, this);
		}
		return tds;
	},

	getTimeStamp: function(){
		return Math.floor(new Date().getTime() * 0.00000001);
	},

	_tempStoreTableData: function(type){
		// caching or clearing table data, depending on the arg
		//
		if(type===true){
			//store indefinitely	
		}else if(type===false){
			// clear object	
			this.tableData = null;
		}else if(type===undefined){
			console.warn("_tempStoreTableData must be passed an argument");	
		}else{ 
			// type is a number/ms
			setTimeout(dojo.hitch(this, function(){
				this.tableData = null;											 
			}), type);
		}
	},

	_tempAvailability: function(type){
			// caching or clearing availability, depending on the arg
		if(type===true){
			//store indefinitely
			this.availableCurrentlySet = true;
		}else if(type===false){
			// clear object	
			this.availableCurrentlySet = false;
		}else if(type===undefined){
			console.warn("_tempAvailability must be passed an argument");	
		}else{ 
			// type is a number/ms
			this.availableCurrentlySet = true;
			setTimeout(dojo.hitch(this, function(){
				this.availableCurrentlySet = false;
			}), type);
		}

	},

	connectTableKeys: function(){
		// summary:
		//		When a table is in focus, start detecting keys
		//		Mainly checking for the TAB key so user can tab 
		//		through a table (blocking the browser's desire to
		//		tab away from teh editor completely)
		if(this.tablesConnected){ return; }
		this.tablesConnected = true;
		var node = (this.editor.iframe) ? this.editor.document : this.editor.editNode;
		this.cnKeyDn = dojo.connect(node, "onkeydown", this, "onKeyDown"); 
		this.cnKeyUp = dojo.connect(node, "onkeyup", this, "onKeyUp");
		this._myListeners.push(dojo.connect(node, "onkeypress", this, "onKeyUp"));
	},

	disconnectTableKeys: function(){
		dojo.disconnect(this.cnKeyDn);
		dojo.disconnect(this.cnKeyUp);
		this.tablesConnected = false;
	},

	onKeyDown: function(evt){
		var key = evt.keyCode;
		if(key == 16){ this.shiftKeyDown = true;}
		if(key == 9) {
			var o = this.getTableInfo();
			// modifying the o.tdIndex in the tableData directly, because we may save it
			// FIXME: tabTo is a global
			o.tdIndex = (this.shiftKeyDown) ? o.tdIndex-1 : tabTo = o.tdIndex+1;
			if(o.tdIndex>=0 && o.tdIndex<o.tds.length){

				this.editor.selectElement(o.tds[o.tdIndex]);

				// we know we are still within a table, so block the need
				//	to run the method
				this.currentlyAvailable = true;
				this._tempAvailability(true);
				// 
				this._tempStoreTableData(true);
				this.stopEvent = true;
			}else{
				//tabbed out of table
				this.stopEvent = false;
				this.onDisplayChanged();
			}
			if(this.stopEvent) {
				dojo.stopEvent(evt);
			}
		}
	},

	onKeyUp: function(evt){
		var key = evt.keyCode;
		if(key == 16){ this.shiftKeyDown = false;}
		if(key == 37 || key == 38 || key == 39 || key == 40 ){
			// user can arrow or tab out of table - need to recheck
			this.onDisplayChanged();
		}
		if(key == 9 && this.stopEvent){ dojo.stopEvent(evt);}
	},

	onDisplayChanged: function(){
		this.currentlyAvailable = false;
		this._tempStoreTableData(false);
		this._tempAvailability(false);
		this.checkAvailable();
	},

	uninitialize: function(editor){
		// summary:
		//		Function to handle cleaning up of connects
		//		and such.  It only finally destroys everything once
		//		all 'references' to it have gone.  As in all plugins
		//		that called init on it destroyed their refs in their 
		//		cleanup calls.
		// editor:
		//		The editor to detach from.
		if(this.editor == editor){
			this.refCount--;
			if(!this.refCount && this.initialized){
				if(this.tablesConnected){
					this.disconnectTableKeys();
				}
				this.initialized = false;
				dojo.forEach(this._myListeners, function(l){
					dojo.disconnect(l);
				});
				delete this._myListeners;
				delete this.editor._tablePluginHandler;
				delete this.editor;
			}
			this.inherited(arguments);
		}
	}
});

dojo.declare("dojox.editor.plugins.TablePlugins",
	dijit._editor._Plugin,
	{
		//summary: 
		//		A collection of Plugins for inserting and modifying tables in the Editor
		//		See end of this document for all avaiable plugs
		//		and dojox/editorPlugins/tests/editorTablePlugs.html for an example
		//
		//		NOT IMPLEMENTED: Not handling cell merge, span or split
		//

		iconClassPrefix: "editorIcon",
		useDefaultCommand: false,
		buttonClass: dijit.form.Button,
		commandName:"",
		label:"",
		alwaysAvailable:false,
		undoEnabled:true,

		onDisplayChanged: function(withinTable){
			// subscribed to from the global object's publish method
			//
			if(!this.alwaysAvailable){
				this.available = withinTable;
				this.button.set('disabled', !this.available);
			}
		},

		setEditor: function(editor){
			this.editor = editor;
			this.editor.customUndo = true;
			this.inherited(arguments);
			this._availableTopic = dojo.subscribe(this.editor.id + "_tablePlugins", this, "onDisplayChanged");
			this.onEditorLoaded();
		},
		onEditorLoaded: function(){
			if(!this.editor._tablePluginHandler){
				// Create it and init it off the editor.  This
				// will create the _tablePluginHandler reference on
				// the dijit.Editor instance.  This avoids a global.
				var tablePluginHandler = new dojox.editor.plugins._TableHandler(); 
				tablePluginHandler.initialize(this.editor);
			}else{
				this.editor._tablePluginHandler.initialize(this.editor);
			}
		},

		selectTable: function(){
			// selects table that is in focus 
			var o = this.getTableInfo();
			if(o && o.tbl){
				dojo.withGlobal(this.editor.window, "selectElement",dijit._editor.selection, [o.tbl]);
			}
		},

		_initButton: function(){
			this.command = this.commandName;

			this.label = this.editor.commands[this.command] = this._makeTitle(this.command);
			this.inherited(arguments);
			delete this.command;

			this.connect(this.button, "onClick", "modTable");

			this.onDisplayChanged(false);
		},

		modTable: function(cmd, args){
			// summary:
			//		Where each plugin performs its action
			//		Note: not using execCommand. In spite of their presence in the 
			//		Editor as query-able plugins, I was not able to find any evidence
			//		that they are supported (especially in NOT IE). If they are 
			//		supported in other browsers, it may help with the undo problem.
			//
			this.begEdit();
			var o = this.getTableInfo();
			var sw = (dojo.isString(cmd))?cmd : this.commandName;
			var r, c, i;
			var adjustColWidth = false;

			if(dojo.isIE){
				// IE can lose selections on focus changes, so focus back
				// in order to restore it.
				this.editor.focus();
			}
			switch(sw){
				case "insertTableRowBefore":
					r = o.tbl.insertRow(o.trIndex);
					for(i=0;i<o.cols;i++){
						c = r.insertCell(-1);
						c.innerHTML = "&nbsp;";
					}
					break;
				case "insertTableRowAfter":
					r = o.tbl.insertRow(o.trIndex+1);
					for(i=0;i<o.cols;i++){
						c = r.insertCell(-1);
						c.innerHTML = "&nbsp;";
					}
					break;
				case "insertTableColumnBefore":
					o.trs.forEach(function(r){
						c = r.insertCell(o.colIndex);
						c.innerHTML = "&nbsp;";
					});
					adjustColWidth = true;
					break;
				case "insertTableColumnAfter":
					o.trs.forEach(function(r){
						c = r.insertCell(o.colIndex+1);
						c.innerHTML = "&nbsp;";
					});
					adjustColWidth = true;
					break;
				case "deleteTableRow":
					o.tbl.deleteRow(o.trIndex);
					break;
				case "deleteTableColumn":
					o.trs.forEach(function(tr){
						tr.deleteCell(o.colIndex);
					});
					adjustColWidth = true;
					break;

				case "modifyTable":
					break;
				case "insertTable":
					break;

			}
			if(adjustColWidth){
				this.makeColumnsEven();
			}
			this.endEdit();
		},

		begEdit: function(){
			if(this.editor._tablePluginHandler.undoEnabled){
				if(this.editor.customUndo){
					this.editor.beginEditing();
				}else{
					this.valBeforeUndo = this.editor.getValue();

				}
			}
		},
		endEdit: function(){
			if(this.editor._tablePluginHandler.undoEnabled){
				if(this.editor.customUndo){
					this.editor.endEditing();
				}else{
					// This code ALMOST works for undo - 
					//	It seems to only work for one step
					//	back in history however
					var afterUndo = this.editor.getValue();
					//this.editor.execCommand("inserthtml", "<p>mike</p>");
					this.editor.setValue(this.valBeforeUndo);
					this.editor.replaceValue(afterUndo);
				}

				this.editor.onDisplayChanged();
			}
		},

		makeColumnsEven: function(){
			//summary:
			//		After changing column amount, change widths to
			//		keep columns even
			//
			// the timeout helps prevent an occasional snafu
			setTimeout(dojo.hitch(this, function(){
				var o = this.getTableInfo(true);
				var w = Math.floor(100/o.cols);
				o.tds.forEach(function(d){
					dojo.attr(d, "width", w+"%");
				});
			}), 10);
		},

		getTableInfo: function(forceNewData){
			// summary:
			//		Gets the table in focus
			//		Collects info on the table - see return params
			//
			return this.editor._tablePluginHandler.getTableInfo(forceNewData);
		},
		_makeTitle: function(str){
			// Parses the commandName into a Title
			//	based on camelCase
			var ns = [];
			dojo.forEach(str, function(c, i){
				if(c.charCodeAt(0)<91 && i>0 && ns[i-1].charCodeAt(0)!=32){
					ns.push(" ");
				}
				if(i===0){ c = c.toUpperCase();}
				ns.push(c);
			});
			return ns.join("");	
		},



		getSelectedCells: function(){
			// summary:
			//		Gets the selected cells from the passed table
			//		Returns: array of TDs or empty array
			var cells = [];
			var tbl = this.getTableInfo().tbl;
			this.editor._tablePluginHandler._prepareTable(tbl);
			var e = this.editor;

			// Lets do this the way IE originally was (Looking up ids).  Walking the selection
			// is inconsistent in the browsers (and painful), so going by ids is simpler.
			var text = dojo.withGlobal(e.window, "getSelectedHtml",dijit._editor.selection, [null]);
			var str = text.match(/id="*\w*"*/g);
			dojo.forEach(str, function(a){
				var id = a.substring(3, a.length);
				if(id.charAt(0) == "\"" && id.charAt(id.length - 1) == "\""){
					id = id.substring(1, id.length - 1);
				}
				var node = e.byId(id);
				if(node && node.tagName.toLowerCase() == "td"){
					cells.push(node);
				}
			}, this);

			if(!cells.length){
				//May just be in a cell (cursor point, or selection in a cell), so look upwards.
				//for a cell container.
				var sel = dijit.range.getSelection(e.window);
				if(sel.rangeCount){
					var r = sel.getRangeAt(0);
					var node = r.startContainer;
					while(node && node != e.editNode && node != e.document){
						if(node.nodeType === 1){
							var tg = node.tagName ? node.tagName.toLowerCase() : "";
							if(tg === "td"){
								return [node];
							}
						}
						node = node.parentNode;
					}
				}
			}
			return cells;
		},

		destroy: function(){
			// summary:
			//		Over-ridden destroy to do some cleanup.
			this.inherited(arguments);
			dojo.unsubscribe(this._availableTopic);

			// Disconnect the editor from the handler
			// to clean up refs.  Moved to using a per-editor
			// 'handler' to avoid collisions on the old global.
			this.editor._tablePluginHandler.uninitialize(this.editor);
		}

	}
);

dojo.declare("dojox.editor.plugins.TableContextMenu",
	dojox.editor.plugins.TablePlugins,
	{
		constructor: function(){
			// summary:
			//		Initialize certain plugins
			//
			this.connect(this, "setEditor", function(editor){
				editor.onLoadDeferred.addCallback(dojo.hitch(this, function() {
					this._createContextMenu();
				}));
				this.button.domNode.style.display = "none";
			});
		},	

		_initButton: function(){
			this.inherited(arguments);
			if(this.commandName=="tableContextMenu"){ this.button.domNode.display = "none";}
		},

		_createContextMenu: function(){
			// summary
			//		Building context menu for right-click shortcuts within a table
			//

			var pMenu = new dijit.Menu({targetNodeIds:[this.editor.iframe]});
			var messages = dojo.i18n.getLocalization("dojox.editor.plugins", "TableDialog", this.lang);
			pMenu.addChild(new dijit.MenuItem({label: messages.selectTableLabel, onClick: dojo.hitch(this, "selectTable")}));
			pMenu.addChild(new dijit.MenuSeparator());

			pMenu.addChild(new dijit.MenuItem({label: messages.insertTableRowBeforeLabel, onClick: dojo.hitch(this, "modTable", "insertTableRowBefore" )}));
			pMenu.addChild(new dijit.MenuItem({label: messages.insertTableRowAfterLabel, onClick: dojo.hitch(this, "modTable", "insertTableRowAfter" )}));
			pMenu.addChild(new dijit.MenuItem({label: messages.insertTableColumnBeforeLabel, onClick: dojo.hitch(this, "modTable", "insertTableColumnBefore" )}));
			pMenu.addChild(new dijit.MenuItem({label: messages.insertTableColumnAfterLabel, onClick: dojo.hitch(this, "modTable", "insertTableColumnAfter" )}));
			pMenu.addChild(new dijit.MenuSeparator());
			pMenu.addChild(new dijit.MenuItem({label: messages.deleteTableRowLabel, onClick: dojo.hitch(this, "modTable", "deleteTableRow" )}));
			pMenu.addChild(new dijit.MenuItem({label: messages.deleteTableColumnLabel, onClick: dojo.hitch(this, "modTable", "deleteTableColumn" )}));

			this.menu = pMenu;
		}		
});

dojo.declare("dojox.editor.plugins.InsertTable",
	dojox.editor.plugins.TablePlugins,
	{
		alwaysAvailable: true,

		modTable: function(){
			var w = new dojox.editor.plugins.EditorTableDialog({});
			w.show();
			var c = dojo.connect(w, "onBuildTable", this, function(obj){
				dojo.disconnect(c);

				var res = this.editor.execCommand('inserthtml', obj.htmlText);

				// commenting this line, due to msg below
				//var td = this.editor.query("td", this.editor.byId(obj.id));

				//HMMMM.... This throws a security error now. didn't used to.
				//this.editor.selectElement(td);
			});
		}	
});

dojo.declare("dojox.editor.plugins.ModifyTable",
	dojox.editor.plugins.TablePlugins,
	{
		modTable: function(){
			if (!this.editor._tablePluginHandler.checkAvailable()) {return;} 
			var o = this.getTableInfo();
			var w = new dojox.editor.plugins.EditorModifyTableDialog({table:o.tbl});
			w.show();
			this.connect(w, "onSetTable", function(color){
				// uhm... not sure whats going on here...
				var o = this.getTableInfo();
				dojo.attr(o.td, "bgcolor", color);
			});
		}	
});

dojo.declare("dojox.editor.plugins.ColorTableCell",
	dojox.editor.plugins.TablePlugins,
	{

		constructor: function(){
			// summary:
			//		Initialize certain plugins
			//
			this.buttonClass = dijit.form.DropDownButton;
			this.dropDown = new dijit.ColorPalette();
			this.connect(this.dropDown, "onChange", function(color){
				this.modTable(null, color);
			});
		},

		_initButton: function(){
			this.command = this.commandName;

			this.label = this.editor.commands[this.command] = this._makeTitle(this.command);
			this.inherited(arguments);
			delete this.command;

			this.onDisplayChanged(false);
		},

		modTable: function(cmd, args){
			// summary
			//	Where each plugin performs its action
			//	Note: not using execCommand. In spite of their presence in the 
			//	Editor as query-able plugins, I was not able to find any evidence
			//	that they are supported (especially in NOT IE). If they are 
			//	supported in other browsers, it may help with the undo problem.
			//
			this.begEdit();
			var o = this.getTableInfo();			
			// The one plugin that really needs use of the very verbose
			//	getSelectedCells()
			var tds = this.getSelectedCells(o.tbl);
			//console.debug("SELECTED CELLS ", tds , " FOR ", o);
			dojo.forEach(tds, function(td){
				dojo.style(td, "backgroundColor", args);				   
			});
			this.endEdit();
		}
});

dojo.provide("dojox.editor.plugins.EditorTableDialog");
dojo.require("dijit.Dialog");
dojo.require("dijit.form.TextBox");
dojo.require("dijit.form.FilteringSelect");
dojo.require("dijit.form.Button");

dojo.declare("dojox.editor.plugins.EditorTableDialog", [dijit.Dialog], {
	// summary:
	//		Dialog box with options for table creation
	//
	baseClass:"EditorTableDialog",

	widgetsInTemplate:true,
	templateString: dojo.cache("dojox.editor.plugins", "resources/insertTable.html"),

	postMixInProperties: function(){
		var messages = dojo.i18n.getLocalization("dojox.editor.plugins", "TableDialog", this.lang);
		dojo.mixin(this, messages);
		this.inherited(arguments);
	},

	postCreate: function(){
		dojo.addClass(this.domNode, this.baseClass); //FIXME - why isn't Dialog accepting the baseClass?
		this.inherited(arguments);
	},

	onInsert: function(){

		var rows =		this.selectRow.get("value") || 1,
			cols =		this.selectCol.get("value") || 1,
			width =		this.selectWidth.get("value"),
			widthType = this.selectWidthType.get("value"),
			border =	this.selectBorder.get("value"),
			pad =		this.selectPad.get("value"),
			space =		this.selectSpace.get("value"),
			_id =		"tbl_"+(new Date().getTime()),
			t = '<table id="'+_id+'"width="'+width+((widthType=="percent")?'%':'')+'" border="'+border+'" cellspacing="'+space+'" cellpadding="'+pad+'">\n';

		for(var r=0;r<rows;r++){
			t += '\t<tr>\n';
			for(var c=0;c<cols;c++){
				t += '\t\t<td width="'+(Math.floor(100/cols))+'%">&nbsp;</td>\n';
			}
			t += '\t</tr>\n';
		}
		t += '</table>';

		this.onBuildTable({htmlText:t, id:_id});
		var cl = dojo.connect(this, "onHide", function(){
			dojo.disconnect(cl);
			var self = this;
			setTimeout(function(){
				self.destroyRecursive();
			}, 10);
		});
		this.hide();
	},

	onCancel: function(){
		// summary:
		//		Function to clean up memory so that the dialog is destroyed 
		//		when closed.
		var c = dojo.connect(this, "onHide", function(){
			dojo.disconnect(c);
			var self = this;
			setTimeout(function(){
				self.destroyRecursive();
			}, 10);
		});
	},

	onBuildTable: function(tableText){
		//stub
	}
});


dojo.provide("dojox.editor.plugins.EditorModifyTableDialog");
dojo.require("dijit.ColorPalette");

dojo.declare("dojox.editor.plugins.EditorModifyTableDialog", [dijit.Dialog], {

	// summary:
	//		Dialog box with options for editing a table
	//

	baseClass:"EditorTableDialog",

	widgetsInTemplate:true,
	table:null, //html table to be modified
	tableAtts:{},
	templateString: dojo.cache("dojox.editor.plugins", "resources/modifyTable.html"),

	postMixInProperties: function(){
		var messages = dojo.i18n.getLocalization("dojox.editor.plugins", "TableDialog", this.lang);
		dojo.mixin(this, messages);
		this.inherited(arguments);
	},

	postCreate: function(){
		dojo.addClass(this.domNode, this.baseClass); //FIXME - why isn't Dialog accepting the baseClass?
		this.inherited(arguments);
		this._cleanupWidgets = [];
		var w1 = new dijit.ColorPalette({});
		this.connect(w1, "onChange", function(color){
			dijit.popup.close(w1);
			this.setBrdColor(color);
		});
		this.connect(w1, "onBlur", function(){
			dijit.popup.close(w1);	
		});
		this.connect(this.borderCol, "click", function(){
			dijit.popup.open({popup:w1, around:this.borderCol});
			w1.focus();
		});
		var w2 = new dijit.ColorPalette({});
		this.connect(w2, "onChange", function(color){
			dijit.popup.close(w2);
			this.setBkColor(color);
		});
		this.connect(w2, "onBlur", function(){
			dijit.popup.close(w2);	
		});
		this.connect(this.backgroundCol, "click", function(){
            dijit.popup.open({popup:w2, around:this.backgroundCol});
			w2.focus();
		});
		this._cleanupWidgets.push(w1);
		this._cleanupWidgets.push(w2);

		this.setBrdColor(dojo.attr(this.table, "bordercolor"));
		this.setBkColor(dojo.attr(this.table, "bgcolor"));
		var w = dojo.attr(this.table, "width");
		if(!w){
			w = this.table.style.width;
		}
		var p = "pixels";
		if(dojo.isString(w) && w.indexOf("%")>-1){
			p = "percent";
			w = w.replace(/%/, "");
		}

		if(w){
			this.selectWidth.set("value", w);
			this.selectWidthType.set("value", p);
		}else{
			this.selectWidth.set("value", "");
			this.selectWidthType.set("value", "percent");
		}

		this.selectBorder.set("value", dojo.attr(this.table, "border"));
		this.selectPad.set("value", dojo.attr(this.table, "cellPadding"));
		this.selectSpace.set("value", dojo.attr(this.table, "cellSpacing"));
		this.selectAlign.set("value", dojo.attr(this.table, "align"));
	},

	setBrdColor: function(color){
		this.brdColor = color;
		dojo.style(this.borderCol, "backgroundColor", color);
	},

	setBkColor: function(color){
		this.bkColor = color;
		dojo.style(this.backgroundCol, "backgroundColor", color);
	},
	onSet: function(){
		dojo.attr(this.table, "borderColor", this.brdColor);
		dojo.attr(this.table, "bgColor", this.bkColor);
		if(this.selectWidth.get("value")){
			// Just in case, remove it from style since we're setting it as a table attribute.
			dojo.style(this.table, "width", "");
			dojo.attr(this.table, "width", (this.selectWidth.get("value") + ((this.selectWidthType.get("value")=="pixels")?"":"%") ));
		}
		dojo.attr(this.table, "border", this.selectBorder.get("value"));
		dojo.attr(this.table, "cellPadding", this.selectPad.get("value"));
		dojo.attr(this.table, "cellSpacing", this.selectSpace.get("value"));
		dojo.attr(this.table, "align", this.selectAlign.get("value"));
		var c = dojo.connect(this, "onHide", function(){
			dojo.disconnect(c);
			var self = this;
			setTimeout(function(){
				self.destroyRecursive();
			}, 10);
		});
		this.hide();
	},

	onCancel: function(){
		// summary:
		//		Function to clean up memory so that the dialog is destroyed 
		//		when closed.
		var c = dojo.connect(this, "onHide", function(){
			dojo.disconnect(c);
			var self = this;
			setTimeout(function(){
				self.destroyRecursive();
			}, 10);
		});
	},

	onSetTable: function(tableText){
		//stub
	},

	destroy: function(){
		// summary:
		//		Cleanup function.
		this.inherited(arguments);
		dojo.forEach(this._cleanupWidgets, function(w){
			if(w && w.destroy){
				w.destroy();
			}
		});
		delete this._cleanupWidgets;
	}
});




dojo.subscribe(dijit._scopeName + ".Editor.getPlugin",null,function(o){
	if(o.plugin){ return; }
	// make first character lower case
	if(o.args && o.args.command){
		var cmd = o.args.command.charAt(0).toLowerCase()+o.args.command.substring(1,o.args.command.length);

		switch(cmd){
			case "insertTableRowBefore":
			case "insertTableRowAfter":
			case "insertTableColumnBefore":
			case "insertTableColumnAfter":
			case "deleteTableRow":
			case "deleteTableColumn":
				o.plugin = new dojox.editor.plugins.TablePlugins({commandName: cmd});
				break;

			case "colorTableCell":
				o.plugin = new dojox.editor.plugins.ColorTableCell({commandName: cmd});
				break;

			case "modifyTable":
				o.plugin = new dojox.editor.plugins.ModifyTable({commandName: cmd});
				break;

			case "insertTable":
				o.plugin = new dojox.editor.plugins.InsertTable({commandName: cmd});
				break;

			case "tableContextMenu":
				o.plugin = new dojox.editor.plugins.TableContextMenu({commandName: cmd});
				break;
		}
	}
});
