"use strict";
(() => {
  var __getOwnPropNames = Object.getOwnPropertyNames;
  var __esm = (fn, res) => function __init() {
    return fn && (res = (0, fn[__getOwnPropNames(fn)[0]])(fn = 0)), res;
  };
  var __commonJS = (cb, mod) => function __require() {
    return mod || (0, cb[__getOwnPropNames(cb)[0]])((mod = { exports: {} }).exports, mod), mod.exports;
  };

  // preview-src/settings.ts
  function getData(key) {
    const element = document.getElementById("vscode-rst-preview-data");
    if (element) {
      const data = element.getAttribute(key);
      if (data) {
        return JSON.parse(data);
      }
    }
    throw new Error(`Could not load data for ${key}`);
  }
  function getSettings() {
    if (cachedSettings) {
      return cachedSettings;
    }
    cachedSettings = getData("data-settings");
    if (cachedSettings) {
      return cachedSettings;
    }
    throw new Error("Could not load settings");
  }
  var cachedSettings;
  var init_settings = __esm({
    "preview-src/settings.ts"() {
      "use strict";
      init_pre();
      cachedSettings = void 0;
    }
  });

  // preview-src/strings.ts
  function getStrings() {
    const store = document.getElementById("vscode-rst-preview-data");
    if (store) {
      const data = store.getAttribute("data-strings");
      if (data) {
        return JSON.parse(data);
      }
    }
    throw new Error("Could not load strings");
  }
  var init_strings = __esm({
    "preview-src/strings.ts"() {
      "use strict";
      init_pre();
    }
  });

  // preview-src/csp.ts
  var CspAlerter;
  var init_csp = __esm({
    "preview-src/csp.ts"() {
      "use strict";
      init_pre();
      init_settings();
      init_strings();
      CspAlerter = class {
        constructor() {
          this.didShow = false;
          this.didHaveCspWarning = false;
          document.addEventListener("securitypolicyviolation", () => {
            this.onCspWarning();
          });
          window.addEventListener("message", (event) => {
            if (event && event.data && event.data.name === "vscode-did-block-svg") {
              this.onCspWarning();
            }
          });
        }
        setPoster(poster) {
          this.messaging = poster;
          if (this.didHaveCspWarning) {
            this.showCspWarning();
          }
        }
        onCspWarning() {
          this.didHaveCspWarning = true;
          this.showCspWarning();
        }
        showCspWarning() {
          const strings = getStrings();
          const settings2 = getSettings();
          if (this.didShow || settings2.disableSecurityWarnings || !this.messaging) {
            return;
          }
          this.didShow = true;
          const notification = document.createElement("a");
          notification.innerText = strings.cspAlertMessageText;
          notification.setAttribute("id", "code-csp-warning");
          notification.setAttribute("title", strings.cspAlertMessageTitle);
          notification.setAttribute("role", "button");
          notification.setAttribute("aria-label", strings.cspAlertMessageLabel);
          notification.onclick = () => {
            this.messaging.postCommand(
              "restructuredtext.showPreviewSecuritySelector",
              [settings2.source]
            );
          };
          document.body.appendChild(notification);
        }
      };
    }
  });

  // preview-src/pre.ts
  var init_pre = __esm({
    "preview-src/pre.ts"() {
      "use strict";
      init_csp();
      window.cspAlerter = new CspAlerter();
    }
  });

  // node_modules/lodash.throttle/index.js
  var require_lodash = __commonJS({
    "node_modules/lodash.throttle/index.js"(exports, module) {
      init_pre();
      var FUNC_ERROR_TEXT = "Expected a function";
      var NAN = 0 / 0;
      var symbolTag = "[object Symbol]";
      var reTrim = /^\s+|\s+$/g;
      var reIsBadHex = /^[-+]0x[0-9a-f]+$/i;
      var reIsBinary = /^0b[01]+$/i;
      var reIsOctal = /^0o[0-7]+$/i;
      var freeParseInt = parseInt;
      var freeGlobal = typeof global == "object" && global && global.Object === Object && global;
      var freeSelf = typeof self == "object" && self && self.Object === Object && self;
      var root = freeGlobal || freeSelf || Function("return this")();
      var objectProto = Object.prototype;
      var objectToString = objectProto.toString;
      var nativeMax = Math.max;
      var nativeMin = Math.min;
      var now = function() {
        return root.Date.now();
      };
      function debounce(func, wait, options) {
        var lastArgs, lastThis, maxWait, result, timerId, lastCallTime, lastInvokeTime = 0, leading = false, maxing = false, trailing = true;
        if (typeof func != "function") {
          throw new TypeError(FUNC_ERROR_TEXT);
        }
        wait = toNumber(wait) || 0;
        if (isObject(options)) {
          leading = !!options.leading;
          maxing = "maxWait" in options;
          maxWait = maxing ? nativeMax(toNumber(options.maxWait) || 0, wait) : maxWait;
          trailing = "trailing" in options ? !!options.trailing : trailing;
        }
        function invokeFunc(time) {
          var args = lastArgs, thisArg = lastThis;
          lastArgs = lastThis = void 0;
          lastInvokeTime = time;
          result = func.apply(thisArg, args);
          return result;
        }
        function leadingEdge(time) {
          lastInvokeTime = time;
          timerId = setTimeout(timerExpired, wait);
          return leading ? invokeFunc(time) : result;
        }
        function remainingWait(time) {
          var timeSinceLastCall = time - lastCallTime, timeSinceLastInvoke = time - lastInvokeTime, result2 = wait - timeSinceLastCall;
          return maxing ? nativeMin(result2, maxWait - timeSinceLastInvoke) : result2;
        }
        function shouldInvoke(time) {
          var timeSinceLastCall = time - lastCallTime, timeSinceLastInvoke = time - lastInvokeTime;
          return lastCallTime === void 0 || timeSinceLastCall >= wait || timeSinceLastCall < 0 || maxing && timeSinceLastInvoke >= maxWait;
        }
        function timerExpired() {
          var time = now();
          if (shouldInvoke(time)) {
            return trailingEdge(time);
          }
          timerId = setTimeout(timerExpired, remainingWait(time));
        }
        function trailingEdge(time) {
          timerId = void 0;
          if (trailing && lastArgs) {
            return invokeFunc(time);
          }
          lastArgs = lastThis = void 0;
          return result;
        }
        function cancel() {
          if (timerId !== void 0) {
            clearTimeout(timerId);
          }
          lastInvokeTime = 0;
          lastArgs = lastCallTime = lastThis = timerId = void 0;
        }
        function flush() {
          return timerId === void 0 ? result : trailingEdge(now());
        }
        function debounced() {
          var time = now(), isInvoking = shouldInvoke(time);
          lastArgs = arguments;
          lastThis = this;
          lastCallTime = time;
          if (isInvoking) {
            if (timerId === void 0) {
              return leadingEdge(lastCallTime);
            }
            if (maxing) {
              timerId = setTimeout(timerExpired, wait);
              return invokeFunc(lastCallTime);
            }
          }
          if (timerId === void 0) {
            timerId = setTimeout(timerExpired, wait);
          }
          return result;
        }
        debounced.cancel = cancel;
        debounced.flush = flush;
        return debounced;
      }
      function throttle2(func, wait, options) {
        var leading = true, trailing = true;
        if (typeof func != "function") {
          throw new TypeError(FUNC_ERROR_TEXT);
        }
        if (isObject(options)) {
          leading = "leading" in options ? !!options.leading : leading;
          trailing = "trailing" in options ? !!options.trailing : trailing;
        }
        return debounce(func, wait, {
          "leading": leading,
          "maxWait": wait,
          "trailing": trailing
        });
      }
      function isObject(value) {
        var type = typeof value;
        return !!value && (type == "object" || type == "function");
      }
      function isObjectLike(value) {
        return !!value && typeof value == "object";
      }
      function isSymbol(value) {
        return typeof value == "symbol" || isObjectLike(value) && objectToString.call(value) == symbolTag;
      }
      function toNumber(value) {
        if (typeof value == "number") {
          return value;
        }
        if (isSymbol(value)) {
          return NAN;
        }
        if (isObject(value)) {
          var other = typeof value.valueOf == "function" ? value.valueOf() : value;
          value = isObject(other) ? other + "" : other;
        }
        if (typeof value != "string") {
          return value === 0 ? value : +value;
        }
        value = value.replace(reTrim, "");
        var isBinary = reIsBinary.test(value);
        return isBinary || reIsOctal.test(value) ? freeParseInt(value.slice(2), isBinary ? 2 : 8) : reIsBadHex.test(value) ? NAN : +value;
      }
      module.exports = throttle2;
    }
  });

  // node_modules/morphdom/dist/morphdom.js
  var require_morphdom = __commonJS({
    "node_modules/morphdom/dist/morphdom.js"(exports, module) {
      "use strict";
      init_pre();
      var DOCUMENT_FRAGMENT_NODE = 11;
      function morphAttrs(fromNode, toNode) {
        var toNodeAttrs = toNode.attributes;
        var attr;
        var attrName;
        var attrNamespaceURI;
        var attrValue;
        var fromValue;
        if (toNode.nodeType === DOCUMENT_FRAGMENT_NODE || fromNode.nodeType === DOCUMENT_FRAGMENT_NODE) {
          return;
        }
        for (var i = toNodeAttrs.length - 1; i >= 0; i--) {
          attr = toNodeAttrs[i];
          attrName = attr.name;
          attrNamespaceURI = attr.namespaceURI;
          attrValue = attr.value;
          if (attrNamespaceURI) {
            attrName = attr.localName || attrName;
            fromValue = fromNode.getAttributeNS(attrNamespaceURI, attrName);
            if (fromValue !== attrValue) {
              if (attr.prefix === "xmlns") {
                attrName = attr.name;
              }
              fromNode.setAttributeNS(attrNamespaceURI, attrName, attrValue);
            }
          } else {
            fromValue = fromNode.getAttribute(attrName);
            if (fromValue !== attrValue) {
              fromNode.setAttribute(attrName, attrValue);
            }
          }
        }
        var fromNodeAttrs = fromNode.attributes;
        for (var d = fromNodeAttrs.length - 1; d >= 0; d--) {
          attr = fromNodeAttrs[d];
          attrName = attr.name;
          attrNamespaceURI = attr.namespaceURI;
          if (attrNamespaceURI) {
            attrName = attr.localName || attrName;
            if (!toNode.hasAttributeNS(attrNamespaceURI, attrName)) {
              fromNode.removeAttributeNS(attrNamespaceURI, attrName);
            }
          } else {
            if (!toNode.hasAttribute(attrName)) {
              fromNode.removeAttribute(attrName);
            }
          }
        }
      }
      var range;
      var NS_XHTML = "http://www.w3.org/1999/xhtml";
      var doc = typeof document === "undefined" ? void 0 : document;
      var HAS_TEMPLATE_SUPPORT = !!doc && "content" in doc.createElement("template");
      var HAS_RANGE_SUPPORT = !!doc && doc.createRange && "createContextualFragment" in doc.createRange();
      function createFragmentFromTemplate(str) {
        var template = doc.createElement("template");
        template.innerHTML = str;
        return template.content.childNodes[0];
      }
      function createFragmentFromRange(str) {
        if (!range) {
          range = doc.createRange();
          range.selectNode(doc.body);
        }
        var fragment = range.createContextualFragment(str);
        return fragment.childNodes[0];
      }
      function createFragmentFromWrap(str) {
        var fragment = doc.createElement("body");
        fragment.innerHTML = str;
        return fragment.childNodes[0];
      }
      function toElement(str) {
        str = str.trim();
        if (HAS_TEMPLATE_SUPPORT) {
          return createFragmentFromTemplate(str);
        } else if (HAS_RANGE_SUPPORT) {
          return createFragmentFromRange(str);
        }
        return createFragmentFromWrap(str);
      }
      function compareNodeNames(fromEl, toEl) {
        var fromNodeName = fromEl.nodeName;
        var toNodeName = toEl.nodeName;
        var fromCodeStart, toCodeStart;
        if (fromNodeName === toNodeName) {
          return true;
        }
        fromCodeStart = fromNodeName.charCodeAt(0);
        toCodeStart = toNodeName.charCodeAt(0);
        if (fromCodeStart <= 90 && toCodeStart >= 97) {
          return fromNodeName === toNodeName.toUpperCase();
        } else if (toCodeStart <= 90 && fromCodeStart >= 97) {
          return toNodeName === fromNodeName.toUpperCase();
        } else {
          return false;
        }
      }
      function createElementNS(name, namespaceURI) {
        return !namespaceURI || namespaceURI === NS_XHTML ? doc.createElement(name) : doc.createElementNS(namespaceURI, name);
      }
      function moveChildren(fromEl, toEl) {
        var curChild = fromEl.firstChild;
        while (curChild) {
          var nextChild = curChild.nextSibling;
          toEl.appendChild(curChild);
          curChild = nextChild;
        }
        return toEl;
      }
      function syncBooleanAttrProp(fromEl, toEl, name) {
        if (fromEl[name] !== toEl[name]) {
          fromEl[name] = toEl[name];
          if (fromEl[name]) {
            fromEl.setAttribute(name, "");
          } else {
            fromEl.removeAttribute(name);
          }
        }
      }
      var specialElHandlers = {
        OPTION: function(fromEl, toEl) {
          var parentNode = fromEl.parentNode;
          if (parentNode) {
            var parentName = parentNode.nodeName.toUpperCase();
            if (parentName === "OPTGROUP") {
              parentNode = parentNode.parentNode;
              parentName = parentNode && parentNode.nodeName.toUpperCase();
            }
            if (parentName === "SELECT" && !parentNode.hasAttribute("multiple")) {
              if (fromEl.hasAttribute("selected") && !toEl.selected) {
                fromEl.setAttribute("selected", "selected");
                fromEl.removeAttribute("selected");
              }
              parentNode.selectedIndex = -1;
            }
          }
          syncBooleanAttrProp(fromEl, toEl, "selected");
        },
        INPUT: function(fromEl, toEl) {
          syncBooleanAttrProp(fromEl, toEl, "checked");
          syncBooleanAttrProp(fromEl, toEl, "disabled");
          if (fromEl.value !== toEl.value) {
            fromEl.value = toEl.value;
          }
          if (!toEl.hasAttribute("value")) {
            fromEl.removeAttribute("value");
          }
        },
        TEXTAREA: function(fromEl, toEl) {
          var newValue = toEl.value;
          if (fromEl.value !== newValue) {
            fromEl.value = newValue;
          }
          var firstChild = fromEl.firstChild;
          if (firstChild) {
            var oldValue = firstChild.nodeValue;
            if (oldValue == newValue || !newValue && oldValue == fromEl.placeholder) {
              return;
            }
            firstChild.nodeValue = newValue;
          }
        },
        SELECT: function(fromEl, toEl) {
          if (!toEl.hasAttribute("multiple")) {
            var selectedIndex = -1;
            var i = 0;
            var curChild = fromEl.firstChild;
            var optgroup;
            var nodeName;
            while (curChild) {
              nodeName = curChild.nodeName && curChild.nodeName.toUpperCase();
              if (nodeName === "OPTGROUP") {
                optgroup = curChild;
                curChild = optgroup.firstChild;
              } else {
                if (nodeName === "OPTION") {
                  if (curChild.hasAttribute("selected")) {
                    selectedIndex = i;
                    break;
                  }
                  i++;
                }
                curChild = curChild.nextSibling;
                if (!curChild && optgroup) {
                  curChild = optgroup.nextSibling;
                  optgroup = null;
                }
              }
            }
            fromEl.selectedIndex = selectedIndex;
          }
        }
      };
      var ELEMENT_NODE = 1;
      var DOCUMENT_FRAGMENT_NODE$1 = 11;
      var TEXT_NODE = 3;
      var COMMENT_NODE = 8;
      function noop() {
      }
      function defaultGetNodeKey(node) {
        if (node) {
          return node.getAttribute && node.getAttribute("id") || node.id;
        }
      }
      function morphdomFactory(morphAttrs2) {
        return function morphdom3(fromNode, toNode, options) {
          if (!options) {
            options = {};
          }
          if (typeof toNode === "string") {
            if (fromNode.nodeName === "#document" || fromNode.nodeName === "HTML" || fromNode.nodeName === "BODY") {
              var toNodeHtml = toNode;
              toNode = doc.createElement("html");
              toNode.innerHTML = toNodeHtml;
            } else {
              toNode = toElement(toNode);
            }
          } else if (toNode.nodeType === DOCUMENT_FRAGMENT_NODE$1) {
            toNode = toNode.firstElementChild;
          }
          var getNodeKey = options.getNodeKey || defaultGetNodeKey;
          var onBeforeNodeAdded = options.onBeforeNodeAdded || noop;
          var onNodeAdded = options.onNodeAdded || noop;
          var onBeforeElUpdated = options.onBeforeElUpdated || noop;
          var onElUpdated = options.onElUpdated || noop;
          var onBeforeNodeDiscarded = options.onBeforeNodeDiscarded || noop;
          var onNodeDiscarded = options.onNodeDiscarded || noop;
          var onBeforeElChildrenUpdated = options.onBeforeElChildrenUpdated || noop;
          var skipFromChildren = options.skipFromChildren || noop;
          var addChild = options.addChild || function(parent, child) {
            return parent.appendChild(child);
          };
          var childrenOnly = options.childrenOnly === true;
          var fromNodesLookup = /* @__PURE__ */ Object.create(null);
          var keyedRemovalList = [];
          function addKeyedRemoval(key) {
            keyedRemovalList.push(key);
          }
          function walkDiscardedChildNodes(node, skipKeyedNodes) {
            if (node.nodeType === ELEMENT_NODE) {
              var curChild = node.firstChild;
              while (curChild) {
                var key = void 0;
                if (skipKeyedNodes && (key = getNodeKey(curChild))) {
                  addKeyedRemoval(key);
                } else {
                  onNodeDiscarded(curChild);
                  if (curChild.firstChild) {
                    walkDiscardedChildNodes(curChild, skipKeyedNodes);
                  }
                }
                curChild = curChild.nextSibling;
              }
            }
          }
          function removeNode(node, parentNode, skipKeyedNodes) {
            if (onBeforeNodeDiscarded(node) === false) {
              return;
            }
            if (parentNode) {
              parentNode.removeChild(node);
            }
            onNodeDiscarded(node);
            walkDiscardedChildNodes(node, skipKeyedNodes);
          }
          function indexTree(node) {
            if (node.nodeType === ELEMENT_NODE || node.nodeType === DOCUMENT_FRAGMENT_NODE$1) {
              var curChild = node.firstChild;
              while (curChild) {
                var key = getNodeKey(curChild);
                if (key) {
                  fromNodesLookup[key] = curChild;
                }
                indexTree(curChild);
                curChild = curChild.nextSibling;
              }
            }
          }
          indexTree(fromNode);
          function handleNodeAdded(el) {
            onNodeAdded(el);
            var curChild = el.firstChild;
            while (curChild) {
              var nextSibling = curChild.nextSibling;
              var key = getNodeKey(curChild);
              if (key) {
                var unmatchedFromEl = fromNodesLookup[key];
                if (unmatchedFromEl && compareNodeNames(curChild, unmatchedFromEl)) {
                  curChild.parentNode.replaceChild(unmatchedFromEl, curChild);
                  morphEl(unmatchedFromEl, curChild);
                } else {
                  handleNodeAdded(curChild);
                }
              } else {
                handleNodeAdded(curChild);
              }
              curChild = nextSibling;
            }
          }
          function cleanupFromEl(fromEl, curFromNodeChild, curFromNodeKey) {
            while (curFromNodeChild) {
              var fromNextSibling = curFromNodeChild.nextSibling;
              if (curFromNodeKey = getNodeKey(curFromNodeChild)) {
                addKeyedRemoval(curFromNodeKey);
              } else {
                removeNode(curFromNodeChild, fromEl, true);
              }
              curFromNodeChild = fromNextSibling;
            }
          }
          function morphEl(fromEl, toEl, childrenOnly2) {
            var toElKey = getNodeKey(toEl);
            if (toElKey) {
              delete fromNodesLookup[toElKey];
            }
            if (!childrenOnly2) {
              if (onBeforeElUpdated(fromEl, toEl) === false) {
                return;
              }
              morphAttrs2(fromEl, toEl);
              onElUpdated(fromEl);
              if (onBeforeElChildrenUpdated(fromEl, toEl) === false) {
                return;
              }
            }
            if (fromEl.nodeName !== "TEXTAREA") {
              morphChildren(fromEl, toEl);
            } else {
              specialElHandlers.TEXTAREA(fromEl, toEl);
            }
          }
          function morphChildren(fromEl, toEl) {
            var skipFrom = skipFromChildren(fromEl);
            var curToNodeChild = toEl.firstChild;
            var curFromNodeChild = fromEl.firstChild;
            var curToNodeKey;
            var curFromNodeKey;
            var fromNextSibling;
            var toNextSibling;
            var matchingFromEl;
            outer:
              while (curToNodeChild) {
                toNextSibling = curToNodeChild.nextSibling;
                curToNodeKey = getNodeKey(curToNodeChild);
                while (!skipFrom && curFromNodeChild) {
                  fromNextSibling = curFromNodeChild.nextSibling;
                  if (curToNodeChild.isSameNode && curToNodeChild.isSameNode(curFromNodeChild)) {
                    curToNodeChild = toNextSibling;
                    curFromNodeChild = fromNextSibling;
                    continue outer;
                  }
                  curFromNodeKey = getNodeKey(curFromNodeChild);
                  var curFromNodeType = curFromNodeChild.nodeType;
                  var isCompatible = void 0;
                  if (curFromNodeType === curToNodeChild.nodeType) {
                    if (curFromNodeType === ELEMENT_NODE) {
                      if (curToNodeKey) {
                        if (curToNodeKey !== curFromNodeKey) {
                          if (matchingFromEl = fromNodesLookup[curToNodeKey]) {
                            if (fromNextSibling === matchingFromEl) {
                              isCompatible = false;
                            } else {
                              fromEl.insertBefore(matchingFromEl, curFromNodeChild);
                              if (curFromNodeKey) {
                                addKeyedRemoval(curFromNodeKey);
                              } else {
                                removeNode(curFromNodeChild, fromEl, true);
                              }
                              curFromNodeChild = matchingFromEl;
                            }
                          } else {
                            isCompatible = false;
                          }
                        }
                      } else if (curFromNodeKey) {
                        isCompatible = false;
                      }
                      isCompatible = isCompatible !== false && compareNodeNames(curFromNodeChild, curToNodeChild);
                      if (isCompatible) {
                        morphEl(curFromNodeChild, curToNodeChild);
                      }
                    } else if (curFromNodeType === TEXT_NODE || curFromNodeType == COMMENT_NODE) {
                      isCompatible = true;
                      if (curFromNodeChild.nodeValue !== curToNodeChild.nodeValue) {
                        curFromNodeChild.nodeValue = curToNodeChild.nodeValue;
                      }
                    }
                  }
                  if (isCompatible) {
                    curToNodeChild = toNextSibling;
                    curFromNodeChild = fromNextSibling;
                    continue outer;
                  }
                  if (curFromNodeKey) {
                    addKeyedRemoval(curFromNodeKey);
                  } else {
                    removeNode(curFromNodeChild, fromEl, true);
                  }
                  curFromNodeChild = fromNextSibling;
                }
                if (curToNodeKey && (matchingFromEl = fromNodesLookup[curToNodeKey]) && compareNodeNames(matchingFromEl, curToNodeChild)) {
                  if (!skipFrom) {
                    addChild(fromEl, matchingFromEl);
                  }
                  morphEl(matchingFromEl, curToNodeChild);
                } else {
                  var onBeforeNodeAddedResult = onBeforeNodeAdded(curToNodeChild);
                  if (onBeforeNodeAddedResult !== false) {
                    if (onBeforeNodeAddedResult) {
                      curToNodeChild = onBeforeNodeAddedResult;
                    }
                    if (curToNodeChild.actualize) {
                      curToNodeChild = curToNodeChild.actualize(fromEl.ownerDocument || doc);
                    }
                    addChild(fromEl, curToNodeChild);
                    handleNodeAdded(curToNodeChild);
                  }
                }
                curToNodeChild = toNextSibling;
                curFromNodeChild = fromNextSibling;
              }
            cleanupFromEl(fromEl, curFromNodeChild, curFromNodeKey);
            var specialElHandler = specialElHandlers[fromEl.nodeName];
            if (specialElHandler) {
              specialElHandler(fromEl, toEl);
            }
          }
          var morphedNode = fromNode;
          var morphedNodeType = morphedNode.nodeType;
          var toNodeType = toNode.nodeType;
          if (!childrenOnly) {
            if (morphedNodeType === ELEMENT_NODE) {
              if (toNodeType === ELEMENT_NODE) {
                if (!compareNodeNames(fromNode, toNode)) {
                  onNodeDiscarded(fromNode);
                  morphedNode = moveChildren(fromNode, createElementNS(toNode.nodeName, toNode.namespaceURI));
                }
              } else {
                morphedNode = toNode;
              }
            } else if (morphedNodeType === TEXT_NODE || morphedNodeType === COMMENT_NODE) {
              if (toNodeType === morphedNodeType) {
                if (morphedNode.nodeValue !== toNode.nodeValue) {
                  morphedNode.nodeValue = toNode.nodeValue;
                }
                return morphedNode;
              } else {
                morphedNode = toNode;
              }
            }
          }
          if (morphedNode === toNode) {
            onNodeDiscarded(fromNode);
          } else {
            if (toNode.isSameNode && toNode.isSameNode(morphedNode)) {
              return;
            }
            morphEl(morphedNode, toNode, childrenOnly);
            if (keyedRemovalList) {
              for (var i = 0, len = keyedRemovalList.length; i < len; i++) {
                var elToRemove = fromNodesLookup[keyedRemovalList[i]];
                if (elToRemove) {
                  removeNode(elToRemove, elToRemove.parentNode, false);
                }
              }
            }
          }
          if (!childrenOnly && morphedNode !== fromNode && fromNode.parentNode) {
            if (morphedNode.actualize) {
              morphedNode = morphedNode.actualize(fromNode.ownerDocument || doc);
            }
            fromNode.parentNode.replaceChild(morphedNode, fromNode);
          }
          return morphedNode;
        };
      }
      var morphdom2 = morphdomFactory(morphAttrs);
      module.exports = morphdom2;
    }
  });

  // preview-src/index.ts
  init_pre();

  // preview-src/activeLineMarker.ts
  init_pre();

  // preview-src/scroll-sync.ts
  init_pre();
  init_settings();
  var codeLineClass = "linemarker";
  var lineNumberPrefix = "linemarker-";
  function clamp(min, max, value) {
    return Math.min(max, Math.max(min, value));
  }
  function clampLine(line) {
    return clamp(0, getSettings().lineCount - 1, line);
  }
  var getCodeLineElements = (() => {
    let elements;
    return () => {
      if (!elements) {
        elements = [{ element: document.body, line: 0 }];
        for (const element of document.getElementsByClassName(
          codeLineClass
        )) {
          for (const className of element.classList) {
            if (className.startsWith(lineNumberPrefix)) {
              const line = +className.substring(
                lineNumberPrefix.length
              );
              if (isNaN(line)) {
                continue;
              }
              if (element.tagName === "CODE" && element.parentElement && element.parentElement.tagName === "PRE") {
                elements.push({
                  element: element.parentElement,
                  line
                });
              } else {
                elements.push({
                  element,
                  line
                });
              }
            }
          }
        }
      }
      return elements;
    };
  })();
  function getElementsForSourceLine(targetLine) {
    const lineNumber = Math.floor(targetLine);
    const lines = getCodeLineElements();
    let previous = lines[0] || null;
    for (const entry of lines) {
      if (entry.line === lineNumber) {
        return { previous: entry, next: void 0 };
      } else if (entry.line > lineNumber) {
        return { previous, next: entry };
      }
      previous = entry;
    }
    return { previous };
  }
  function getLineElementsAtPageOffset(offset) {
    const lines = getCodeLineElements();
    const position = offset - window.scrollY;
    let lo = -1;
    let hi = lines.length - 1;
    while (lo + 1 < hi) {
      const mid = Math.floor((lo + hi) / 2);
      const bounds = getElementBounds(lines[mid]);
      if (bounds.top + bounds.height >= position) {
        hi = mid;
      } else {
        lo = mid;
      }
    }
    const hiElement = lines[hi];
    const hiBounds = getElementBounds(hiElement);
    if (hi >= 1 && hiBounds.top > position) {
      const loElement = lines[lo];
      return { previous: loElement, next: hiElement };
    }
    if (hi > 1 && hi < lines.length && hiBounds.top + hiBounds.height > position) {
      return { previous: hiElement, next: lines[hi + 1] };
    }
    return { previous: hiElement };
  }
  function getElementBounds({ element }) {
    const myBounds = element.getBoundingClientRect();
    const codeLineChild = element.querySelector(`.${codeLineClass}`);
    if (codeLineChild) {
      const childBounds = codeLineChild.getBoundingClientRect();
      const height = Math.max(1, childBounds.top - myBounds.top);
      return {
        top: myBounds.top,
        height
      };
    }
    return myBounds;
  }
  function scrollToRevealSourceLine(line) {
    if (!getSettings().scrollPreviewWithEditor) {
      return;
    }
    if (line <= 0) {
      window.scroll(window.scrollX, 0);
      return;
    }
    const { previous, next } = getElementsForSourceLine(line);
    if (!previous) {
      return;
    }
    let scrollTo = 0;
    const rect = getElementBounds(previous);
    const previousTop = rect.top;
    if (next && next.line !== previous.line) {
      const betweenProgress = (line - previous.line) / (next.line - previous.line);
      const elementOffset = next.element.getBoundingClientRect().top - previousTop;
      scrollTo = previousTop + betweenProgress * elementOffset;
    } else {
      const progressInElement = line - Math.floor(line);
      scrollTo = previousTop + rect.height * progressInElement;
    }
    window.scroll(window.scrollX, Math.max(1, window.scrollY + scrollTo));
  }
  function getEditorLineNumberForPageOffset(offset) {
    const { previous, next } = getLineElementsAtPageOffset(offset);
    if (previous) {
      const previousBounds = getElementBounds(previous);
      const offsetFromPrevious = offset - window.scrollY - previousBounds.top;
      if (next) {
        const progressBetweenElements = offsetFromPrevious / (getElementBounds(next).top - previousBounds.top);
        const line = previous.line + progressBetweenElements * (next.line - previous.line);
        return clampLine(line);
      } else {
        const progressWithinElement = offsetFromPrevious / previousBounds.height;
        const line = previous.line + progressWithinElement;
        return clampLine(line);
      }
    }
    return null;
  }

  // preview-src/activeLineMarker.ts
  var ActiveLineMarker = class {
    onDidChangeTextEditorSelection(line) {
      const { previous } = getElementsForSourceLine(line);
      this._update(previous && previous.element);
    }
    _update(before) {
      this._unmarkActiveElement(this._current);
      this._markActiveElement(before);
      this._current = before;
    }
    _unmarkActiveElement(element) {
      if (!element) {
        return;
      }
      element.className = element.className.replace(
        /\bcode-active-line\b/g,
        ""
      );
    }
    _markActiveElement(element) {
      if (!element) {
        return;
      }
      element.className += " code-active-line";
    }
  };

  // preview-src/events.ts
  init_pre();
  function onceDocumentLoaded(f) {
    if (document.readyState === "loading") {
      document.addEventListener("DOMContentLoaded", f);
    } else {
      f();
    }
  }

  // preview-src/messaging.ts
  init_pre();
  init_settings();
  var createPosterForVsCode = (vscode2) => {
    return new class {
      postMessage(type, body) {
        vscode2.postMessage({
          type,
          source: getSettings().source,
          body
        });
      }
      postCommand(command, args) {
        this.postMessage("command", { command, args });
      }
    }();
  };

  // preview-src/index.ts
  init_settings();
  var throttle = require_lodash();
  var morphdom = require_morphdom();
  var scrollDisabled = true;
  var marker = new ActiveLineMarker();
  var settings = getSettings();
  var vscode = acquireVsCodeApi();
  var state = getData("data-state");
  vscode.setState(state);
  var messaging = createPosterForVsCode(vscode);
  window.cspAlerter.setPoster(messaging);
  onceDocumentLoaded(() => {
    if (settings.scrollPreviewWithEditor) {
      setTimeout(() => {
        const initialLine = +settings.line;
        if (!isNaN(initialLine)) {
          scrollDisabled = true;
          scrollToRevealSourceLine(initialLine);
        }
      }, 0);
    }
  });
  var onUpdateView = (() => {
    const doScroll = throttle((line) => {
      scrollDisabled = true;
      scrollToRevealSourceLine(line);
    }, 50);
    return (line, settings2) => {
      if (!isNaN(line)) {
        settings2.line = line;
        doScroll(line);
      }
    };
  })();
  window.addEventListener(
    "resize",
    () => {
      scrollDisabled = true;
    },
    true
  );
  window.addEventListener(
    "message",
    (event) => {
      if (event.data.source !== settings.source) {
        return;
      }
      switch (event.data.type) {
        case "onDidChangeTextEditorSelection":
          marker.onDidChangeTextEditorSelection(event.data.line);
          break;
        case "updateContent":
          const parser = new DOMParser();
          const target_doc = parser.parseFromString(
            event.data.content,
            "text/html"
          );
          morphdom(document.body, target_doc.body, {
            childrenOnly: true
          });
          window.document.dispatchEvent(
            new Event("DOMContentLoaded", {
              bubbles: true,
              cancelable: true
            })
          );
          if ("MathJax" in window) {
            MathJax.typeset();
          }
          break;
        case "updateView":
          onUpdateView(event.data.line, settings);
          break;
      }
    },
    false
  );
  document.addEventListener("dblclick", (event) => {
    if (!settings.doubleClickToSwitchToEditor) {
      return;
    }
    for (let node = event.target; node; node = node.parentNode) {
      if (node.tagName === "A") {
        return;
      }
    }
    const offset = event.pageY;
    const line = getEditorLineNumberForPageOffset(offset);
    if (typeof line === "number" && !isNaN(line)) {
      messaging.postMessage("didClick", { line: Math.floor(line) });
    }
  });
  function getLineNumber(elem) {
    const classname = elem.className;
    const line_classname = classname.split(" ").find((cn) => cn.startsWith("linemarker-"));
    const marker_vs = line_classname.split("-")[1];
    return marker_vs;
  }
  document.addEventListener(
    "click",
    (event) => {
      if (!event) {
        return;
      }
      if (event.ctrlKey) {
        const pos = event.clientY + window.scrollY;
        const elems = document.getElementsByClassName("linemarker");
        if (elems.length === 0) {
          return;
        }
        let found = false;
        let marker_start = 0;
        let marker_end = 0;
        for (let i = 0; i < elems.length; i++) {
          const elem = elems.item(i);
          const role = elem?.parentElement?.getAttribute("role");
          if (role === "navigation") {
            continue;
          }
          const pos_elem = elem.getBoundingClientRect().top + window.scrollY;
          if (pos_elem > pos) {
            if (i > 0) {
              marker_start = getLineNumber(elems.item(i - 1));
              marker_end = getLineNumber(elem);
            } else {
              marker_start = -1;
              marker_end = getLineNumber(elem);
            }
            found = true;
            break;
          }
        }
        if (found === false) {
          marker_start = getLineNumber(elems.item(elems.length - 1));
          marker_end = -1;
        }
        vscode.postMessage({
          command: "show_range",
          start: marker_start,
          end: marker_end
        });
        event.stopPropagation();
        event.stopImmediatePropagation();
        event.preventDefault();
        return false;
      }
      let node = event.target;
      while (node) {
        if (node.tagName && node.tagName === "A" && node.href) {
          if (node.getAttribute("href").startsWith("#")) {
            break;
          }
          if (node.href.startsWith("file://") || node.href.startsWith("vscode-resource:")) {
            const [path, fragment] = node.href.replace(/^(file:\/\/|vscode-resource:)/i, "").split("#");
            messaging.postCommand("_rst.openDocumentLink", [
              { path, fragment }
            ]);
            event.preventDefault();
            event.stopPropagation();
            break;
          }
          break;
        }
        node = node.parentNode;
      }
    },
    true
  );
  if (settings.scrollEditorWithPreview) {
    window.addEventListener(
      "scroll",
      throttle(() => {
        if (scrollDisabled) {
          scrollDisabled = false;
        } else {
          const line = getEditorLineNumberForPageOffset(window.scrollY);
          if (typeof line === "number" && !isNaN(line)) {
            messaging.postMessage("revealLine", { line });
          }
        }
      }, 50)
    );
  }
  window.onload = function() {
    vscode.postMessage({ command: "load" });
  };
})();
