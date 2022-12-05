"use strict";
exports.__esModule = true;
var query_any_every_helpers_1 = require("../query_any_every_helpers");
describe('injectValueAtPath', function () {
    it('correctly injects and appends ::any', function () {
        var array = ['modelName', 'secondModel', ['thirdModel', '::any'], '::any'];
        var injected = query_any_every_helpers_1.injectValueAtPath(array, [2, 1], ['new', 'tuple']);
        expect(array).toEqual([
            'modelName',
            'secondModel',
            ['thirdModel', ['new', 'tuple'], '::any'],
            '::any'
        ]);
        expect(injected).toEqual(true);
    });
    it('correctly injects and appends ::every', function () {
        var array = [
            'modelName',
            'secondModel',
            ['thirdModel', '::every'],
            '::every'
        ];
        var injected = query_any_every_helpers_1.injectValueAtPath(array, [2, 1], ['new', 'tuple']);
        expect(array).toEqual([
            'modelName',
            'secondModel',
            ['thirdModel', ['new', 'tuple'], '::every'],
            '::every'
        ]);
        expect(injected).toEqual(true);
    });
    it('unpacks the new value (instead of injecting) if no ::any being replaced', function () {
        var array = ['modelName', ['something', ['nested']]];
        var injected = query_any_every_helpers_1.injectValueAtPath(array, [1, 1, 1], ['new', 'tuple']);
        expect(array).toEqual([
            'modelName',
            ['something', ['nested', 'new', 'tuple']]
        ]);
        expect(injected).toEqual(false);
    });
});
