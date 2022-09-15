from typing import List,Dict
import torch
from collections import OrderedDict

LENTHS={'angstrom':1.8897161646321,
        'bohr':1}

def _check_unit(unit:str,unitdict):
    if unit in unitdict:
        pass
    else:
        raise ValueError(f'invalid unit. allows: {unitdict.keys()}')

def _to_unit(tensor,ori_unit,to_unit,units):
    scale=units[to_unit]/units[ori_unit]
    return tensor/scale

class RawGeo:
    def __init__(self) -> None:
        self.segments:Dict[str,List[List[str]]]=OrderedDict()

    @classmethod 
    def from_in(cls,in_file):
        cls=cls()
        key='holder'
        with open(in_file,'r') as f:
            for i in f.readlines():
                if i.startswith('%'):
                    key=i[1:-1]
                    assert key not in cls.segments.keys(),f'key redefined: {key}'
                    cls.segments[key]=[]
                else:
                    if len(i.strip())>0:
                        cls.segments[key].append(i.split())
        return cls
    @classmethod
    def from_xyz(cls,xyz_file):
        cls=cls()
        with open(xyz_file,'r') as f:
            for i,j in enumerate(f.readlines):
                if i==0:
                    cls.segments['ATOM_NUMBER']=j.split()
                elif i==1:
                    cls.segments['META_INFO']=j.split()
                elif i==2:
                    cls.segments['ATOMIC_POSTION']=[j.split()]
                else:
                    cls.segments['ATOMIC_POSTION'].append(j.split())

    def __getitem__(self,__k)->List[List[str]]:
        return self.segments[__k]

    def out(self,outfile:str):
        if outfile.endswith('xyz'):
            pass
        elif outfile.endswith('in'):
            pass
        else:
            raise NotImplementedError

    def _out_xyz(self,outfile:str):
        pass

    def _out_in(self,outfile:str):
        pass


class Coordination:
    def __init__(self,raw:List[List[str]],unit:str='angstrom') -> None:
        self.U=LENTHS
        _check_unit(unit,self.U)
        cord_str=[i[1:] for i in raw]
        cord_float = [ [ float(_) for _ in i] for i in cord_str]
        self.cord=torch.Tensor(cord_float)
        self.cord=self.cord.type(torch.float32)
        self.atom=[i[0] for i in raw]
        self.unit=unit

    def to(self,unit:str)->None:
        _check_unit(unit,self.U)
        self.cord=_to_unit(self.cord,self.unit,unit,self.U)
        # scale=self.U[unit]/self.U[self.unit]
        # self.cord=self.cord/scale
        self.unit=unit


class CellParameter:
    def __init__(self,raw:List[List[str]],unit:str='angstrom') -> None:
        self.U=LENTHS
        _check_unit(unit,self.U)
        self.unit=unit
        self.base=OrderedDict()
        def raw2base(r):
            bl=[float(i) for i in r]
            bt=torch.Tensor(bl)
            return bt.type(torch.float32)  
        self.base['x']=raw2base(raw[0])
        self.base['y']=raw2base(raw[1])
        self.base['z']=raw2base(raw[2])
        self.fullbase=torch.cat([i.view(1,-1) for i in self.base.values()],dim=0)

    def to(self,unit:str)->None:
        _check_unit(unit)
